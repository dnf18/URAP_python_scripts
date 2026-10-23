import os
import json
import smtplib
from email.message import EmailMessage
from pathlib import Path
from xml.sax.saxutils import escape

from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors

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
        doc = SimpleDocTemplate(str(output_path), pagesize=letter)
        styles = getSampleStyleSheet()
        cell = styles["BodyText"]
        cell.fontSize = 8
        cell.leading = 10
        story  = []

        # Title
        story.append(Paragraph("<b>TGRS MEGAlib Validation Report</b>", styles["Title"]))
        story.append(Spacer(1, 14))

        # Configuration Table (long values wrap inside their cells)
        story.append(Paragraph("<b>Configuration Summary</b>", styles["Heading2"]))
        if self.config:
            rows = [[Paragraph("<b>Parameter</b>", cell), Paragraph("<b>Value</b>", cell)]]
            for k, v in _flatten(self.config):
                rows.append([Paragraph(escape(k), cell), Paragraph(escape(v), cell)])
            config_table = Table(rows, colWidths=[200, 300], repeatRows=1)
            config_table.setStyle(TableStyle([
                ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("GRID", (0, 0), (-1, -1), 0.5, colors.black),
            ]))
            story.append(config_table)
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
        for name, stats in all_stats.items():
            decides = (name == DECIDING_SECTION)
            heading = name if decides else f"{name} (information only)"
            story.append(Paragraph(f"<b>{heading}</b>", styles["Heading2"]))
            story.append(Spacer(1, 6))

            # Plot (the plots are 2:1)
            plot_path = plot_files.get(name)
            if plot_path and plot_path.exists():
                story.append(Image(str(plot_path), width=460, height=230))
                story.append(Spacer(1, 8))
            else:
                story.append(Paragraph("Plot not found.", styles["BodyText"]))

            # Stats table
            stat_data = [["Metric", "Value"]]
            for k, v in stats.items():
                if k == "overlay_plot":
                    continue
                stat_data.append([str(k), str(round(v, 6) if isinstance(v, float) else v)])

            stat_table = Table(stat_data, colWidths=[280, 220])
            stat_table.setStyle(TableStyle([
                ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
                ("GRID", (0, 0), (-1, -1), 0.5, colors.black),
                ("ALIGN", (0, 0), (-1, -1), "LEFT"),
            ]))
            story.append(stat_table)
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
        story.append(Paragraph(
            f"<b>Overall Validation Status: <font color='{overall_color}'>{overall_status}</font></b>",
            styles["Heading2"]
        ))
        story.append(Paragraph("Decided by the A versus B comparison only.", styles["BodyText"]))

        doc.build(story)
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
