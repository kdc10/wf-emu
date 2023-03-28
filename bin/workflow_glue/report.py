"""Create workflow report."""
import json

from ezcharts.components import fastcat
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Tabs
from ezcharts.layout.snippets.table import DataTable
import pandas as pd
import workflow_glue.report_utils.report_utils as report_utils

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    report = labs.LabsReport(
        "Workflow Emu report", "wf-emu",
        args.params, args.versions)

    with open(args.metadata) as metadata:
        sample_details = sorted([
            {
                'sample': d['alias'],
                'type': d['type'],
                'barcode': d['barcode']
            } for d in json.load(metadata)
        ], key=lambda d: d["sample"])

    # Add a section with statistic per sample
    if args.stats:
        with report.add_section("Read summary", "Read summary"):
            fastcat.SeqSummary(args.stats)

    # Add a section with main EMU results
    samples_tables = {table.split('_')[0]: table for table in args.rel_abun}
    logger.info(f"{samples_tables}.")
    with report.add_section("Results", "Results"):
        tabs = Tabs()
        # 2.1. Table
        if len(samples_tables) == 1:
            sample_id = args.rel_abun[0].split('_')[0]
            df = pd.read_csv(args.rel_abun[0], index_col=0, sep='\t').round(3)
            report_utils.DataTable.from_pandas(
                df, export=True, file_name=f'wf-emu-{sample_id}_rel-abundance.tsv')
        else:
            # add drowpdown tabs
            with tabs.add_dropdown_menu('Results', change_header=False):
                for sample_id, rel_table in sorted(samples_tables.items()):
                    with tabs.add_dropdown_tab(sample_id):
                        df = pd.read_csv(rel_table, index_col=0, sep='\t').round(3)
                        report_utils.DataTable.from_pandas(
                            df, export=True, file_name=f'wf-emu-{sample_id}_rel-abundance.tsv')

    with report.add_section("Metadata", "Metadata"):
        tabs = Tabs()
        for d in sample_details:
            with tabs.add_tab(d["sample"]):
                df = pd.DataFrame.from_dict(d, orient="index", columns=["Value"])
                df.index.name = "Key"
                DataTable.from_pandas(df)

    report.write(args.report)
    logger.info(f"Report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("report", help="Report output file")
    parser.add_argument("--stats", nargs='*', help="Fastcat per-read stats file(s).")
    parser.add_argument("--rel_abun", nargs='*', help="Relative abundance.")
    parser.add_argument(
        "--metadata", default='metadata.json',
        help="sample metadata")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A JSON file containing the workflow parameter key/values")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    return parser
