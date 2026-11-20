import os
import json
import smtplib
from email.message import EmailMessage
from pathlib import Path
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.enums import TA_RIGHT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
from reportlab.platypus.tableofcontents import TableOfContents

"""
This is the TGRS Reporter Class

This is where we attach all relevant data/graphs onto a PDF named TGRS_Validation_Report.pdf
It takes the 5 plots and 5 JSON stat files produced by the TGRS comparator and assembles them
into a single PDF report.

The overall PASS/FAIL is decided ONLY by the "Sim A vs Sim B (Combined)" section, i.e. by
the build-to-build comparison. The comparisons with the measurements are shown as
information: they depend on the simulation statistics and on physics that is not part of
a build regression test.
"""

# The one section that decides the overall status
DECIDING_SECTION = "Sim A vs Sim B (Combined)"


class _ReportDocument(SimpleDocTemplate):
    def afterFlowable(self, flowable):
        entry = getattr(flowable, "toc_entry", None)
        if entry is not None:
            title, key = entry
            self.notify("TOCEntry", (0, escape(title), self.page, key))
            self.canv.addOutlineEntry(title, key)


def _section_heading(number, title, style):
    key = f"section-{number}"
    label = f"{number}. {title}"
    heading = Paragraph(f'<a name="{key}"/>{escape(label)}', style)
    heading.toc_entry = (label, key)
    return heading


def _page_number(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 9)
    canvas.setFillColor(colors.HexColor("#475569"))
    canvas.drawRightString(doc.pagesize[0] - doc.rightMargin, 29, str(doc.page))
    canvas.restoreState()


def _report_table(headers, values, column_widths, cell, numeric):
    rows = [[Paragraph(f"<b>{escape(label)}</b>", cell) for label in headers]]
    for label, value in values:
        is_number = isinstance(value, (int, float)) and not isinstance(value, bool)
        display = format(value, ".6g") if isinstance(value, float) else str(value)
        rows.append([
            Paragraph(escape(str(label)), cell),
            Paragraph(escape(display), numeric if is_number else cell),
        ])

    table = Table(rows, colWidths=column_widths, repeatRows=1, hAlign="LEFT")
    table.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f6f8fa")]),
        ("LINEABOVE", (0, 0), (-1, 0), 0.7, colors.HexColor("#64748b")),
        ("LINEBELOW", (0, 0), (-1, 0), 0.5, colors.HexColor("#94a3b8")),
        ("LINEBELOW", (0, -1), (-1, -1), 0.7, colors.HexColor("#64748b")),
    ]))
    return table


def _flatten(data, prefix=""):
    """Turns nested dicts into [["a / b", "value"], ...] rows."""
    rows = []
    for k, v in data.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            rows += _flatten(v, key + " / ")
        else:
            rows.append([key, str(v)])
    return rows


class TGRS_Reporter:

    def __init__(self, config_json: str = None, output_dir: str = "./reports", email_recipients=None):

        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Load run configuration (run_config.json from the comparator, or a steering config)
        self.config = {}
        if config_json and os.path.exists(config_json):
            with open(config_json, "r") as f:
                self.config = json.load(f)

        self.email_recipients = email_recipients or []

    def load_stats(self, stats_dir: str = "./results"):
        """
        Loads all 5 stat JSON files produced by the TGRS comparator.
        Returns a dict of {name: stats_dict}
        """
        stats_dir = Path(stats_dir)
        stat_files = {
            "Sim A vs Measurements":       stats_dir / "reference_stats.json",
            "Sim B vs Measurements":       stats_dir / "test_stats.json",
            "Ratio A vs Measurements":     stats_dir / "reference_ratio_stats.json",
            "Ratio B vs Measurements":     stats_dir / "test_ratio_stats.json",
            DECIDING_SECTION:              stats_dir / "combined_stats.json",
        }

        all_stats = {}
        for name, path in stat_files.items():
            if path.exists():
                with open(path) as f:
                    all_stats[name] = json.load(f)
            else:
                print(f"[Reporter] Stats file not found: {path}")

        return all_stats

    def generate_pdf(self, stats_dir: str = "./results", output_name: str = "TGRS_Validation_Report.pdf"):
        """
        Creates PDF report with all 5 plots and stat tables.
        """

        stats_dir  = Path(stats_dir)
        output_path = self.output_dir / output_name
        doc = _ReportDocument(str(output_path), pagesize=letter,
                              leftMargin=48, rightMargin=48, topMargin=48, bottomMargin=54,
                              title="TGRS MEGAlib Validation Report")
        styles = getSampleStyleSheet()
        section_style = ParagraphStyle("ReportSection", parent=styles["Heading2"], keepWithNext=True)
        cell = ParagraphStyle("TableCell", parent=styles["BodyText"], fontSize=9, leading=12)
        numeric = ParagraphStyle("NumericCell", parent=cell, alignment=TA_RIGHT)
        table_width = doc.width - 12  # Default frame padding on both sides.
        story = []

        # Title
        story.append(Paragraph("<b>TGRS MEGAlib Validation Report</b>", styles["Title"]))
        story.append(Paragraph("Contents", styles["Heading2"]))
        contents = TableOfContents()
        contents.levelStyles = [ParagraphStyle("ContentsEntry", parent=styles["BodyText"],
                                               fontSize=10, leading=13, spaceBefore=0, rightIndent=24)]
        contents.dotsMinLevel = 0
        story.extend([contents, Spacer(1, 12)])

        story.append(_section_heading(1, "Configuration Summary", section_style))
        if self.config:
            story.append(_report_table(("Parameter", "Value"), _flatten(self.config),
                                       (table_width * 0.39, table_width * 0.61), cell, numeric))
        else:
            story.append(Paragraph("No configuration provided.", styles["BodyText"]))
        story.append(Spacer(1, 12))

        # Load all stats
        all_stats = self.load_stats(stats_dir)

        # Plot files
        plot_files = {
            "Sim A vs Measurements":       stats_dir / "plot1_sim_a.png",
            "Sim B vs Measurements":       stats_dir / "plot2_sim_b.png",
            "Ratio A vs Measurements":     stats_dir / "plot3_ratio_a.png",
            "Ratio B vs Measurements":     stats_dir / "plot4_ratio_b.png",
            DECIDING_SECTION:              stats_dir / "plot5_all.png",
        }

        deciding_pass = None   # None = the deciding section is missing

        # One section per comparison
        for number, (name, stats) in enumerate(all_stats.items(), start=2):
            decides = (name == DECIDING_SECTION)
            heading = name if decides else f"{name} (information only)"
            story.append(_section_heading(number, heading, section_style))

            # Plot (the plots are 2:1)
            plot_path = plot_files.get(name)
            if plot_path and plot_path.exists():
                story.append(Image(str(plot_path), width=460, height=230))
                story.append(Spacer(1, 8))
            else:
                story.append(Paragraph("Plot not found.", styles["BodyText"]))

            # Residual plot between the two simulations (deciding section only)
            if decides:
                resid = stats_dir / "plot6_residual_ab.png"
                if resid.exists():
                    story.append(Image(str(resid), width=460, height=268))
                    story.append(Spacer(1, 8))

            # Stats table
            stat_rows = [(key, value) for key, value in stats.items() if key != "overlay_plot"]
            story.append(_report_table(("Metric", "Value"), stat_rows,
                                       (table_width * 0.55, table_width * 0.45), cell, numeric))
            story.append(Spacer(1, 6))

            # Pass/Fail per section
            passed = bool(stats.get("pass", False))
            if decides:
                deciding_pass = passed
            status = "PASS" if passed else "FAIL"
            color  = "green" if passed else "red"
            note = "" if decides else " (information only - does not affect the overall status)"
            story.append(Paragraph(f"Status: <font color='{color}'><b>{status}</b></font>{note}",
                                   styles["BodyText"]))
            story.append(Spacer(1, 16))

        # Overall status: build-to-build comparison only
        if deciding_pass is None:
            overall_status, overall_color = "UNKNOWN (A vs B comparison missing)", "orange"
        else:
            overall_status = "PASS" if deciding_pass else "FAIL"
            overall_color = "green" if deciding_pass else "red"
        story.append(_section_heading(len(all_stats) + 2, "Overall Validation Status", section_style))
        story.append(Paragraph(
            f"<font color='{overall_color}'><b>{overall_status}</b></font>. "
            "Decided by the A versus B comparison only.", styles["BodyText"]
        ))

        doc.multiBuild(story, onFirstPage=_page_number, onLaterPages=_page_number)
        print(f"PDF report generated at: {output_path}")
        return output_path

    def send_email(self, subject: str, body: str, attachment_path: str,
                   sender_email: str, sender_password: str):
        """
        Emails the PDF report to all recipients.
        """

        if not self.email_recipients:
            print("No email recipients specified. Skipping email step.")
            return

        msg = EmailMessage()
        msg["Subject"] = subject
        msg["From"]    = sender_email
        msg["To"]      = ", ".join(self.email_recipients)
        msg.set_content(body)

        with open(attachment_path, "rb") as f:
            msg.add_attachment(
                f.read(),
                maintype="application",
                subtype="pdf",
                filename=os.path.basename(attachment_path)
            )

        try:
            with smtplib.SMTP_SSL("smtp.gmail.com", 465) as smtp:
                smtp.login(sender_email, sender_password)
                smtp.send_message(msg)
            print("Email sent successfully to:", ", ".join(self.email_recipients))
        except Exception as e:
            print("Failed to send email:", e)


# --- Main ---
if __name__ == "__main__":
    reporter = TGRS_Reporter(
        config_json="results/run_config.json",
        output_dir="./reports"
    )

    pdf_path = reporter.generate_pdf(stats_dir="./results")

    # Optional email
    # reporter.send_email(
    #     subject="TGRS Validation Report",
    #     body="Please find the TGRS validation report attached.",
    #     attachment_path=str(pdf_path),
    #     sender_email="you@gmail.com",
    #     sender_password="yourpassword"
    # )
