#!/usr/bin/env python3
"""
EnteroScope Main Orchestrator – Temporary Directory Version
All module writes occur in a temporary directory; final results are copied to the user output.
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School – Department of Medical Biochemistry
Version: 1.1.0
Date: 2026
"""

import os
import sys
import glob
import argparse
import subprocess
import shutil
import tempfile
import logging
import traceback
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional

__version__ = "1.1.0"

# =============================================================================
# ANSI Color Definitions
# =============================================================================
class Color:
    """ANSI escape codes for coloured terminal output."""
    RESET = '\033[0m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    BRIGHT_BLACK = '\033[90m'
    BRIGHT_RED = '\033[91m'
    BRIGHT_GREEN = '\033[92m'
    BRIGHT_YELLOW = '\033[93m'
    BRIGHT_BLUE = '\033[94m'
    BRIGHT_MAGENTA = '\033[95m'
    BRIGHT_CYAN = '\033[96m'
    BRIGHT_WHITE = '\033[97m'


# =============================================================================
# Banner and Message Display
# =============================================================================
class EnteroScopeBanner:
    """Handles banners, quotes, and coloured console messages."""

    def __init__(self):
        self.quotes = [
            {"text": "Enterobacter cloacae complex is a leading cause of hospital-acquired infections, including pneumonia, UTIs, and bloodstream infections.", "author": "CDC"},
            {"text": "Carbapenem-resistant Enterobacteriaceae (CRE) are 'nightmare bacteria' with mortality rates up to 50%.", "author": "CDC"},
            {"text": "Enterobacter species are intrinsically resistant to ampicillin and amoxicillin due to constitutive AmpC beta-lactamase.", "author": "Unknown"},
            {"text": "The OXA-48 carbapenemase is increasingly reported in Enterobacter cloacae, often carried on epidemic plasmids.", "author": "Unknown"},
            {"text": "Biofilm formation by Enterobacter cloacae on medical devices facilitates persistent infections and resistance spread.", "author": "Unknown"},
            {"text": "Enterobacter cloacae can acquire colistin resistance via mcr genes or chromosomal mutations in pmrAB.", "author": "Brown Beckley"},
            {"text": "Tigecycline resistance in Enterobacter is often mediated by tet(X) variants or efflux pumps like AdeABC.", "author": "Unknown"},
            {"text": "Efflux pumps (AcrAB-TolC, OqxAB, AdeFGH) contribute to multidrug resistance in Enterobacter cloacae.", "author": "Unknown"},
            {"text": "The IncL/M plasmid carrying blaOXA-48 has become a global vehicle for carbapenem resistance in Enterobacter cloacae.", "author": "Unknown"},
            {"text": "MLST of Enterobacter cloacae reveals population structure, with ST78, ST114, and ST171 being common outbreak clones.", "author": "Unknown"},
            {"text": "Science is organised knowledge.", "author": "Herbert Spencer"},
            {"text": "The science of today is the technology of tomorrow.", "author": "Edward Teller"},
            {"text": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie"},
            {"text": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson"},
            {"text": "Genomics is a lens on biology.", "author": "Eric Lander"},
            {"text": "Sequence today, understand tomorrow.", "author": "Anonymous"},
            {"text": "Microbes rule the world.", "author": "Paul Stamets"},
            {"text": "In every drop, a universe.", "author": "Antonie van Leeuwenhoek"},
            {"text": "Evolution in a petri dish.", "author": "Richard Lenski"},
            {"text": "Science knows no country, because knowledge belongs to humanity.", "author": "Louis Pasteur"},
            {"text": "Imagination is more important than knowledge.", "author": "Albert Einstein"},
            {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein"},
            {"text": "We cannot solve our problems with the same thinking we used when we created them.", "author": "Albert Einstein"},
            {"text": "Life is like riding a bicycle. To keep your balance, you must keep moving.", "author": "Albert Einstein"},
            {"text": "The only source of knowledge is experience.", "author": "Albert Einstein"},
            {"text": "The greatest danger to our future is apathy.", "author": "Jane Goodall"},
            {"text": "What you do makes a difference, and you have to decide what kind of difference you want to make.", "author": "Jane Goodall"},
            {"text": "The best way to predict the future is to create it.", "author": "Peter Drucker"},
            {"text": "In the middle of difficulty lies opportunity.", "author": "Albert Einstein"},
            {"text": "The beautiful thing about learning is that nobody can take it away from you.", "author": "B.B. King"},
            {"text": "Science is a way of thinking much more than it is a body of knowledge.", "author": "Carl Sagan"},
            {"text": "We are a way for the cosmos to know itself.", "author": "Carl Sagan"},
            {"text": "The universe is not required to be in perfect harmony with human ambition.", "author": "Carl Sagan"},
            {"text": "Somewhere, something incredible is waiting to be known.", "author": "Carl Sagan"},
            {"text": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson"},
            {"text": "The universe is under no obligation to make sense to you.", "author": "Neil deGrasse Tyson"},
            {"text": "We are all connected; To each other, biologically. To the earth, chemically. To the rest of the universe, atomically.", "author": "Neil deGrasse Tyson"},
            {"text": "The nitrogen in our DNA, the calcium in our teeth, the iron in our blood, the carbon in our apple pies were made in the interiors of collapsing stars.", "author": "Carl Sagan"},
            {"text": "The most exciting phrase to hear in science, the one that heralds new discoveries, is not 'Eureka!' but 'That's funny...'", "author": "Isaac Asimov"},
            {"text": "The saddest aspect of life right now is that science gathers knowledge faster than society gathers wisdom.", "author": "Isaac Asimov"},
            {"text": "Any sufficiently advanced technology is indistinguishable from magic.", "author": "Arthur C. Clarke"},
            {"text": "The only way of discovering the limits of the possible is to venture a little way past them into the impossible.", "author": "Arthur C. Clarke"},
            {"text": "There are two great forces in this world: the force of nature and the human spirit.", "author": "TEDx talk"},
            {"text": "Science is a way of life. It is not a collection of facts but a way of thinking.", "author": "Albert Einstein"},
            {"text": "The most beautiful thing we can experience is the mysterious. It is the source of all true art and science.", "author": "Albert Einstein"},
            {"text": "I have no special talent. I am only passionately curious.", "author": "Albert Einstein"},
            {"text": "The scientists of today think deeply instead of clearly. One must be sane to think clearly, but one can think deeply and be quite insane.", "author": "Nikola Tesla"},
            {"text": "If you want to find the secrets of the universe, think in terms of energy, frequency and vibration.", "author": "Nikola Tesla"},
            {"text": "The present is theirs; the future, for which I really worked, is mine.", "author": "Nikola Tesla"},
        ]

    def display_banner(self, show_quote: bool = True, show_author: bool = True) -> None:
        """Print the EnteroScope ASCII banner."""
        banner = f"""
{Color.BOLD}{Color.BRIGHT_CYAN}
███████╗███╗   ██╗████████╗███████╗██████╗  ██████╗ ███████╗ ██████╗  ██████╗ ██████╗ ███████╗
██╔════╝████╗  ██║╚══██╔══╝██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔════╝ ██╔═══██╗██╔══██╗██╔════╝
█████╗  ██╔██╗ ██║   ██║   █████╗  ██████╔╝██║   ██║███████╗██║      ██║   ██║██████╔╝█████╗  
██╔══╝  ██║╚██╗██║   ██║   ██╔══╝  ██╔══██╗██║   ██║╚════██║██║      ██║   ██║██╔═══╝ ██╔══╝  
███████╗██║ ╚████║   ██║   ███████╗██║  ██║╚██████╔╝███████║╚██████╗ ╚██████╔╝██║     ███████╗
╚══════╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝ ╚═════╝  ╚═════╝ ╚═╝     ╚══════╝
{Color.RESET}
{Color.BRIGHT_YELLOW}EnteroScope v{__version__} – Enterobacter cloacae Complex Genomic Analysis Pipeline{Color.RESET}
{Color.CYAN}Author: Brown Beckley | GitHub: https://github.com/bbeckley-hub/enteroscope{Color.RESET}
{Color.CYAN}Affiliation: University of Ghana Medical School – Department of Medical Biochemistry{Color.RESET}
"""
        print(banner)
        if show_quote:
            self.display_random_quote(show_author)

    def display_random_quote(self, show_author: bool = True) -> None:
        """Print a random scientific quote."""
        import random
        if not self.quotes:
            return
        q = random.choice(self.quotes)
        print(f"{Color.DIM}{Color.WHITE}{'─' * 80}{Color.RESET}")
        print(f"{Color.BRIGHT_CYAN}💡 Scientific Insight:{Color.RESET}")
        print(f"{Color.BRIGHT_WHITE}   \"{q['text']}\"{Color.RESET}")
        if show_author and 'author' in q:
            print(f"{Color.BRIGHT_YELLOW}   — {q['author']}{Color.RESET}")
        print(f"{Color.DIM}{Color.WHITE}{'─' * 80}{Color.RESET}")
        print()

    def display_startup_sequence(self) -> None:
        """Display the initial banner and a quote."""
        self.display_banner(show_quote=True, show_author=True)

    def display_info(self, msg: str) -> None:
        """Print an informational message."""
        print(f"{Color.BRIGHT_CYAN}[INFO]{Color.RESET} {msg}")

    def display_success(self, msg: str) -> None:
        """Print a success message."""
        print(f"{Color.BRIGHT_GREEN}✓{Color.RESET} {msg}")

    def display_warning(self, msg: str) -> None:
        """Print a warning message."""
        print(f"{Color.BRIGHT_YELLOW}⚠️{Color.RESET} {msg}")

    def display_error(self, msg: str) -> None:
        """Print an error message."""
        print(f"{Color.BRIGHT_RED}✗{Color.RESET} {msg}")

    def display_module_header(self, title: str, subtitle: str = "") -> None:
        """Print a section header for a module."""
        print()
        print(f"{Color.BOLD}{Color.BRIGHT_BLUE}{'=' * 80}{Color.RESET}")
        print(f"{Color.BOLD}{Color.BRIGHT_CYAN}{' ' * 20}{title}{Color.RESET}")
        if subtitle:
            print(f"{Color.DIM}{Color.WHITE}{' ' * 22}{subtitle}{Color.RESET}")
        print(f"{Color.BOLD}{Color.BRIGHT_BLUE}{'=' * 80}{Color.RESET}")
        print()

    def display_footer(self, analysis_time: str, samples_processed: int) -> None:
        """Print the final summary footer."""
        print(f"\n{Color.BOLD}{Color.BRIGHT_BLUE}{'='*80}{Color.RESET}")
        print(f"{Color.BOLD}{Color.BRIGHT_GREEN}🎉 Analysis Complete!{Color.RESET}")
        print(f"{Color.CYAN}   Time elapsed: {analysis_time}{Color.RESET}")
        print(f"{Color.CYAN}   Samples processed: {samples_processed}{Color.RESET}")
        print(f"{Color.BOLD}{Color.BRIGHT_BLUE}{'='*80}{Color.RESET}")
        print()


# =============================================================================
# Coloured Help Formatter for argparse
# =============================================================================
class ColoredHelpFormatter(argparse.HelpFormatter):
    """Custom argparse help formatter with ANSI colours."""

    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    CYAN = '\033[96m'
    BOLD = '\033[1m'
    RESET = '\033[0m'

    def _format_usage(self, usage, actions, groups, prefix):
        usage_str = super()._format_usage(usage, actions, groups, prefix)
        if usage_str:
            usage_str = self.GREEN + self.BOLD + "usage: " + self.RESET + usage_str
        return usage_str

    def _format_action(self, action):
        action_str = super()._format_action(action)
        if not action_str:
            return action_str
        lines = action_str.split('\n')
        colored_lines = []
        for line in lines:
            if line.strip():
                if line.lstrip().startswith('-'):
                    parts = line.split('  ', 1)
                    if len(parts) == 2:
                        options = parts[0].strip()
                        help_text = parts[1]
                        colored_line = f"  {self.CYAN}{options}{self.RESET}  {help_text}"
                    else:
                        colored_line = f"  {self.CYAN}{line.strip()}{self.RESET}"
                else:
                    colored_line = f"  {self.YELLOW}{line}{self.RESET}"
                colored_lines.append(colored_line)
            else:
                colored_lines.append(line)
        return '\n'.join(colored_lines)

    def _format_text(self, text):
        if not text:
            return text
        return f"{self.BLUE}{text}{self.RESET}"

    def start_section(self, heading):
        heading = f"{self.BOLD}{self.GREEN}{heading}{self.RESET}"
        super().start_section(heading)


# =============================================================================
# Main Orchestrator Class
# =============================================================================
class EnteroScopeOrchestrator:
    """
    Orchestrates the EnteroScope analysis pipeline.

    All modules are run inside temporary directories. Input FASTA files are copied
    there, and results are copied back to the final output directory only upon
    successful completion. This approach keeps the working environment clean and
    is suitable for HPC and Dockerised deployments.
    """

    def __init__(self):
        self.banner = EnteroScopeBanner()
        self.base_dir = Path(__file__).parent
        self.user_output_dir = None
        self.logger = None
        self.keep_temp = False

    def setup_logging(self, output_dir: Path) -> None:
        """Configure logging to a file in the output directory."""
        log_file = output_dir / "enteroscope_run.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[logging.FileHandler(log_file, mode='w')]
        )
        self.logger = logging.getLogger("EnteroScope")
        self.logger.info(f"Logging to {log_file}")
        self.user_output_dir = output_dir

    def get_module_path(self, module_name: str) -> Path:
        """
        Return the path to a module directory.

        The function first checks the installed share directory (if any) and
        otherwise falls back to the local 'modules' subdirectory.
        """
        if hasattr(sys, 'prefix'):
            share_path = Path(sys.prefix) / "share" / "enteroscope" / "modules" / module_name
            if share_path.exists():
                return share_path
        return self.base_dir / "modules" / module_name

    def run_module_in_temp(self, module_name: str, fasta_files: List[Path],
                           cmd_str: str, result_subdir: Optional[str] = None,
                           extra_files: Optional[List[str]] = None,
                           extra_files_from_subdir: Optional[List[str]] = None) -> bool:
        """
        Run a module inside a temporary directory.

        Parameters
        ----------
        module_name : str
            Name of the module folder (e.g., "qc_module").
        fasta_files : List[Path]
            List of FASTA files to copy into the temporary workspace.
        cmd_str : str
            Shell command to execute inside the temporary directory.
        result_subdir : str, optional
            Subdirectory within the temporary directory containing results to copy.
        extra_files : List[str], optional
            Additional file names (relative to the temporary directory) to copy
            to the final output directory.
        extra_files_from_subdir : List[str], optional
            Additional file names (relative to result_subdir) to copy
            to the final output directory after result_subdir is copied.

        Returns
        -------
        bool
            True if the module ran successfully and results were copied, False otherwise.
        """
        module_orig = self.get_module_path(module_name)
        if not module_orig.exists():
            self.banner.display_error(f"Module directory not found: {module_orig}")
            return False

        temp_dir = tempfile.mkdtemp(prefix=f"enteroscope_{module_name}_")
        self.logger.info(f"Temporary directory for {module_name}: {temp_dir}")

        try:
            # Copy module contents and FASTA files
            shutil.copytree(module_orig, Path(temp_dir) / module_name, dirs_exist_ok=True)
            for f in fasta_files:
                shutil.copy2(f, Path(temp_dir) / f.name)

            pattern = self.get_file_pattern(fasta_files)
            self.banner.display_info(f"Running {module_name} analysis with pattern: {pattern}")

            # Execute the command
            result = subprocess.run(cmd_str, shell=True, cwd=temp_dir, capture_output=True, text=True)
            if result.stdout:
                self.logger.info(f"STDOUT from {module_name}:\n{result.stdout}")
            if result.stderr:
                self.logger.warning(f"STDERR from {module_name}:\n{result.stderr}")

            if result.returncode != 0:
                self.logger.error(f"{module_name} failed with return code {result.returncode}")
                return False

            # Copy results
            if result_subdir:
                src = Path(temp_dir) / result_subdir
                if src.exists():
                    dst = self.user_output_dir / result_subdir
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst)
                    self.logger.info(f"Results copied to {dst}")
                    # Copy extra files from the result subdirectory
                    if extra_files_from_subdir:
                        for fname in extra_files_from_subdir:
                            file_src = src / fname
                            if file_src.exists():
                                shutil.copy2(file_src, self.user_output_dir / fname)
                                self.logger.info(f"Copied {fname} from result subdir to output directory")
                            else:
                                self.logger.warning(f"File not found in result subdir: {file_src}")

            if extra_files:
                for fname in extra_files:
                    src = Path(temp_dir) / fname
                    if src.exists():
                        shutil.copy2(src, self.user_output_dir / fname)
                        self.logger.info(f"Copied {fname} to output directory")

            return True

        except Exception as e:
            self.logger.error(f"Exception in {module_name}: {e}\n{traceback.format_exc()}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    # --------------------------------------------------------------------------
    # Module‑specific runners
    # --------------------------------------------------------------------------
    def run_fasta_qc_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run the FASTA QC module."""
        pattern = self.get_file_pattern(fasta_files)
        cmd = f"{sys.executable} qc_module/enteroscope_fasta_qc.py {pattern} -o fasta_qc_results -c {threads}"
        return self.run_module_in_temp("qc_module", fasta_files, cmd, "fasta_qc_results")

    def run_mlst_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """
        Run the MLST module exactly as the manual command:
        cd into mlst_module, run mlst_module.py with -db db -sc bin.
        Input pattern is adjusted to point to the parent directory (where FASTA files are).
        Output is written to ../mlst_results (which is the temp root's mlst_results).
        """
        pattern_unquoted = self.get_file_pattern(fasta_files).strip('"')
        # The FASTA files are in the temp root, which is the parent of mlst_module.
        # So we use "../" + pattern_unquoted to match them from inside mlst_module.
        parent_pattern = f"../{pattern_unquoted}"
        cmd = (
            f"cd mlst_module && "
            f"{sys.executable} mlst_module.py "
            f"-i \"{parent_pattern}\" "
            f"-o ../mlst_results "
            f"-db db "
            f"-sc bin "
            f"--batch"
        )
        return self.run_module_in_temp("mlst_module", fasta_files, cmd, "mlst_results")

    def run_abricate_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int,
                              min_identity: Optional[float] = None,
                              min_coverage: Optional[float] = None) -> bool:
        """
        Run the ABRicate module.

        Parameters
        ----------
        min_identity : float, optional
            Minimum identity percentage for ABRicate.
        min_coverage : float, optional
            Minimum coverage percentage for ABRicate.
        """
        pattern = self.get_file_pattern(fasta_files)
        cmd = f"{sys.executable} abricate_module/enteroscope_abricate.py {pattern} --cpus {threads}"
        if min_identity is not None:
            cmd += f" --min-identity {min_identity}"
        if min_coverage is not None:
            cmd += f" --min-coverage {min_coverage}"
        return self.run_module_in_temp("abricate_module", fasta_files, cmd, "enteroscope_abricate_results")

    def run_amr_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int,
                         min_identity: Optional[float] = None,
                         min_coverage: Optional[float] = None,
                         skip_mutations: bool = False,
                         force_update: bool = False) -> bool:
        """
        Run the AMRfinderPlus module.

        Parameters
        ----------
        min_identity : float, optional
            Minimum identity (0..1) for AMR hits.
        min_coverage : float, optional
            Minimum coverage (0..1) for AMR hits.
        skip_mutations : bool
            If True, disable point mutation reporting.
        force_update : bool
            If True, force an AMR database update before analysis.
        """
        if not self.ensure_amr_database():
            self.banner.display_error("AMR database is missing and could not be updated automatically.")
            return False

        pattern = self.get_file_pattern(fasta_files)
        cmd = f"{sys.executable} amr_module/enteroscope_amrfinder.py {pattern} --cpus {threads}"
        if min_identity is not None:
            cmd += f" --min-identity {min_identity}"
        if min_coverage is not None:
            cmd += f" --min-coverage {min_coverage}"
        if skip_mutations:
            cmd += " --skip-mutations"
        if force_update:
            self.banner.display_info("Forcing AMR database update before analysis...")
            self.update_amr_database(force=True)

        # The mutation files are generated inside the result_subdir (enteroscope_amrfinder_results)
        extra_files_from_subdir = ["mutation_summary.html", "mutation_summary.tsv", "mutation_master_summary.json"]
        return self.run_module_in_temp("amr_module", fasta_files, cmd, "enteroscope_amrfinder_results",
                                       extra_files_from_subdir=extra_files_from_subdir)

    def run_summary_analysis(self, output_dir: Path) -> bool:
        """
        Run the Ultimate Reporter module.

        Required HTML files are collected from the output directory, copied to a
        temporary workspace, and the reporter is executed there. The final reports
        are then copied back.
        """
        temp_dir = tempfile.mkdtemp(prefix="enteroscope_summary_")
        temp_path = Path(temp_dir)
        self.logger.info(f"Temporary directory for summary reports: {temp_dir}")

        try:
            summary_module_orig = self.get_module_path("summary_module")
            shutil.copytree(summary_module_orig, temp_path, dirs_exist_ok=True)

            # List of (subdirectory, filename) pairs to copy
            required_htmls = [
                ("mlst_results", "mlst_summary.html"),
                ("fasta_qc_results", "FASTA_QC_summary.html"),
                ("enteroscope_amrfinder_results", "enteroscope_amrfinder_summary_report.html"),
                ("enteroscope_amrfinder_results", "mutation_summary.html"),
                ("enteroscope_abricate_results", "enteroscope_card_summary_report.html"),
                ("enteroscope_abricate_results", "enteroscope_ncbi_summary_report.html"),
                ("enteroscope_abricate_results", "enteroscope_resfinder_summary_report.html"),
                ("enteroscope_abricate_results", "enteroscope_vfdb_summary_report.html"),
                ("enteroscope_abricate_results", "enteroscope_argannot_summary_report.html"),
                ("enteroscope_abricate_results", "enteroscope_megares_summary_report.html"),
                ("enteroscope_abricate_results", "enteroscope_ecoli_vf_summary_report.html"),
                ("enteroscope_abricate_results", "enteroscope_bacmet2_summary_report.html"),
                ("enteroscope_abricate_results", "enteroscope_plasmidfinder_summary_report.html"),
                ("enteroscope_abricate_results", "enteroscope_ecoh_summary_report.html"),
            ]
            for subdir, filename in required_htmls:
                src = output_dir / subdir / filename
                if src.exists():
                    shutil.copy2(src, temp_path / filename)
                    self.logger.info(f"Copied {filename} to temporary summary_module directory")
                else:
                    self.logger.warning(f"Required HTML not found: {src}")

            # Copy mutation summary files - try multiple locations
            mutation_files = ["mutation_summary.html", "mutation_summary.tsv", "mutation_master_summary.json"]
            sources_to_check = [
                output_dir,  # top-level output directory
                output_dir / "enteroscope_amrfinder_results",  # AMR results subdirectory
            ]
            for fname in mutation_files:
                copied = False
                for src_dir in sources_to_check:
                    src = src_dir / fname
                    if src.exists():
                        shutil.copy2(src, temp_path / fname)
                        self.logger.info(f"Copied {fname} from {src_dir} to temporary summary_module directory")
                        copied = True
                        break
                if not copied:
                    self.logger.warning(f"Mutation file {fname} not found in any location (will be ignored)")

            # Also copy any additional summary files from the root output
            for f in output_dir.glob("*.html"):
                if f.name not in [h[1] for h in required_htmls] and f.name not in mutation_files:
                    shutil.copy2(f, temp_path / f.name)
                    self.logger.info(f"Copied additional HTML: {f.name}")

            self.banner.display_info("Running ultimate reporter...")
            cmd = [sys.executable, "enteroscope_ultimate_reporter.py", "-i", "."]
            self.logger.info(f"Running command: {' '.join(cmd)} (cwd={temp_path})")
            result = subprocess.run(cmd, cwd=temp_path, capture_output=True, text=True)

            # Log stdout and stderr
            if result.stdout:
                self.logger.info(f"STDOUT from ultimate reporter:\n{result.stdout}")
            if result.stderr:
                self.logger.warning(f"STDERR from ultimate reporter:\n{result.stderr}")

            if result.returncode != 0:
                self.logger.error(f"Ultimate reporter failed with return code {result.returncode}")
                self.banner.display_error(f"Ultimate reporter failed with exit code {result.returncode}")
                return False

            self.banner.display_success("Ultimate reporter completed successfully!")
            src_dir = temp_path / "ENTERO_ULTIMATE_REPORTS"
            if src_dir.exists():
                dst_dir = self.user_output_dir / "ENTERO_ULTIMATE_REPORTS"
                if dst_dir.exists():
                    shutil.rmtree(dst_dir)
                shutil.copytree(src_dir, dst_dir)
                self.logger.info(f"Ultimate reports copied to {dst_dir}")
                return True
            else:
                self.banner.display_warning("Ultimate reports directory not found.")
                return False

        except Exception as e:
            self.logger.error(f"Summary reports exception: {e}\n{traceback.format_exc()}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_path, ignore_errors=True)

    # --------------------------------------------------------------------------
    # AMR database handling
    # --------------------------------------------------------------------------
    def update_amr_database(self, force: bool = False) -> bool:
        """Update the AMRfinderPlus database using the AMR module."""
        amr_module_path = self.get_module_path("amr_module")
        amr_script = amr_module_path / "enteroscope_amrfinder.py"
        if not amr_script.exists():
            self.banner.display_error(f"AMR script not found at: {amr_script}")
            return False
        self.banner.display_info("Updating AMRfinderPlus database...")
        flag = "--force-update" if force else "--update-db"
        cmd = [sys.executable, str(amr_script), flag]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=amr_module_path)
        if result.returncode == 0:
            self.banner.display_success("AMR database updated successfully.")
            # Get version but don't print path to terminal – only log
            version_cmd = [sys.executable, str(amr_script), "--db-version"]
            version_result = subprocess.run(version_cmd, capture_output=True, text=True, cwd=amr_module_path)
            if version_result.returncode == 0:
                version_info = version_result.stdout.strip()
                self.logger.info(f"New database version: {version_info}")
                # Only show version, not the path
                self.banner.display_info(f"Database version: {version_info}")
            return True
        else:
            self.banner.display_error("AMR database update failed.")
            if result.stderr:
                self.logger.error(result.stderr)
            return False

    def ensure_amr_database(self) -> bool:
        """Check if the AMR database exists; if not, attempt to update it."""
        amr_module_path = self.get_module_path("amr_module")
        amr_script = amr_module_path / "enteroscope_amrfinder.py"
        if not amr_script.exists():
            self.banner.display_error("AMR script not found, cannot check database.")
            return False
        cmd = [sys.executable, str(amr_script), "--db-version"]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=amr_module_path)
        if result.returncode == 0 and "Unknown" not in result.stdout and "No database" not in result.stdout:
            # Extract version from stdout (first line) and show only version, not path
            version_line = result.stdout.strip().split('\n')[0]
            if "Database version:" in version_line:
                version_info = version_line.replace("Database version:", "").strip()
            else:
                version_info = version_line
            self.banner.display_success(f"AMR database ready (version: {version_info})")
            self.logger.info(f"AMR database present: {result.stdout.strip()}")
            return True
        else:
            self.banner.display_warning("AMR database not found or outdated. Attempting automatic update...")
            return self.update_amr_database(force=False)

    # --------------------------------------------------------------------------
    # Helper methods
    # --------------------------------------------------------------------------
    def find_fasta_files(self, input_path: str) -> List[Path]:
        """
        Locate FASTA files based on the input pattern.

        Supports glob patterns, single files, or directories.
        """
        self.banner.display_info(f"Searching for files with pattern: {input_path}")
        if '*' in input_path or '?' in input_path:
            matched_files = glob.glob(input_path)
            fasta_files = [Path(f) for f in matched_files if Path(f).is_file() and
                           f.lower().endswith(('.fna', '.fasta', '.fa', '.fn')) and
                           not Path(f).name.startswith('.')]
            self.banner.display_success(f"Found {len(fasta_files)} FASTA files")
            return sorted(fasta_files)
        input_path_obj = Path(input_path)
        if input_path_obj.is_file() and input_path_obj.suffix.lower() in ['.fna', '.fasta', '.fa', '.fn']:
            self.banner.display_success(f"Found single FASTA file: {input_path_obj.name}")
            return [input_path_obj]
        if input_path_obj.is_dir():
            patterns = [f"{input_path}/*.fna", f"{input_path}/*.fasta", f"{input_path}/*.fa", f"{input_path}/*.fn"]
            fasta_files = []
            for pattern in patterns:
                for file_path in glob.glob(pattern):
                    path = Path(file_path)
                    if path.is_file() and not path.name.startswith('.'):
                        fasta_files.append(path)
            fasta_files = sorted(list(set(fasta_files)))
            if fasta_files:
                self.banner.display_success(f"Found {len(fasta_files)} FASTA files in directory")
            else:
                self.banner.display_warning(f"No FASTA files found in directory: {input_path}")
            return fasta_files
        self.banner.display_error(f"Input path not found: {input_path}")
        return []

    def get_file_pattern(self, fasta_files: List[Path]) -> str:
        """Return a glob pattern string representing the file extensions of the given files."""
        if not fasta_files:
            return '"*.fna"'
        extensions = set(f.suffix.lower() for f in fasta_files)
        if len(extensions) == 1:
            ext = list(extensions)[0]
            return f'"*{ext}"'
        return '"*"'

    # --------------------------------------------------------------------------
    # Main pipeline driver
    # --------------------------------------------------------------------------
    def run_complete_analysis(
        self,
        input_path: str,
        output_dir: str,
        threads: int = 2,
        skip_modules: Optional[Dict[str, bool]] = None,
        skip_summary: bool = False,
        update_amr_db_only: bool = False,
        amr_min_identity: Optional[float] = None,
        amr_min_coverage: Optional[float] = None,
        amr_skip_mutations: bool = False,
        amr_force_update: bool = False,
        abricate_min_identity: Optional[float] = None,
        abricate_min_coverage: Optional[float] = None,
        clean_output: bool = False
    ) -> None:
        """
        Execute the full EnteroScope analysis pipeline.

        Parameters
        ----------
        input_path : str
            Path or glob pattern to input FASTA files.
        output_dir : str
            Destination directory for all results.
        threads : int, default 2
            Number of CPU threads for modules.
        skip_modules : dict, optional
            Keys: 'qc', 'mlst', 'abricate', 'amr'. Values True to skip.
        skip_summary : bool
            If True, skip the Ultimate Reporter.
        update_amr_db_only : bool
            If True, only update the AMR database and exit.
        amr_min_identity : float, optional
            Minimum identity (0..1) for AMR hits.
        amr_min_coverage : float, optional
            Minimum coverage (0..1) for AMR hits.
        amr_skip_mutations : bool
            If True, disable point mutation reporting in AMR.
        amr_force_update : bool
            If True, force AMR database update before AMR analysis.
        abricate_min_identity : float, optional
            Minimum identity percentage for ABRicate.
        abricate_min_coverage : float, optional
            Minimum coverage percentage for ABRicate.
        clean_output : bool
            If True, delete the output directory before starting.
        """
        if skip_modules is None:
            skip_modules = {}
        if update_amr_db_only:
            self.update_amr_database(force=False)
            return

        start_time = datetime.now()
        try:
            self.banner.display_startup_sequence()
            output_path = Path(output_dir)
            if clean_output and output_path.exists():
                shutil.rmtree(output_path)
            output_path.mkdir(parents=True, exist_ok=True)
            self.setup_logging(output_path)

            fasta_files = self.find_fasta_files(input_path)
            if not fasta_files:
                self.banner.display_error("No FASTA files found! Analysis stopped.")
                return
            extensions = set(f.suffix.lower() for f in fasta_files)
            self.banner.display_success(f"Starting analysis of {len(fasta_files)} samples")
            self.banner.display_info(f"File formats detected: {', '.join(extensions)}")

            # Create output subdirectories (will be filled by copy steps)
            for subdir in ["fasta_qc_results", "mlst_results", "enteroscope_abricate_results",
                           "enteroscope_amrfinder_results", "ENTERO_ULTIMATE_REPORTS"]:
                (output_path / subdir).mkdir(exist_ok=True)

            # Show analysis plan
            self.banner.display_module_header("Analysis Plan", "Modules to be executed")
            analyses = [
                ("FASTA QC", not skip_modules.get('qc', False)),
                ("MLST", not skip_modules.get('mlst', False)),
                ("ABRicate", not skip_modules.get('abricate', False)),
                ("AMR", not skip_modules.get('amr', False)),
                ("Ultimate Reporter", not skip_summary),
            ]
            for name, enabled in analyses:
                status = "✅ ENABLED" if enabled else "⏸️  SKIPPED"
                print(f"   {status} - {name}")
            sys.stdout.flush()

            # ================================================================
            # Sequential execution: QC, then MLST
            # ================================================================
            if not skip_modules.get('qc', False):
                self.banner.display_module_header("FASTA QC", "Sequence Quality Control & Assembly Statistics")
                qc_success = self.run_fasta_qc_analysis(fasta_files, output_path, threads)
                if qc_success:
                    self.banner.display_success("✅ FASTA QC completed")
                else:
                    self.banner.display_error("❌ FASTA QC failed")
                self.banner.display_random_quote()
                print()

            if not skip_modules.get('mlst', False):
                self.banner.display_module_header("MLST", "Multi‑Locus Sequence Typing")
                mlst_success = self.run_mlst_analysis(fasta_files, output_path, threads)
                if mlst_success:
                    self.banner.display_success("✅ MLST completed")
                else:
                    self.banner.display_error("❌ MLST failed")
                self.banner.display_random_quote()
                print()

            # ================================================================
            # Sequential second batch: ABRicate, then AMR
            # ================================================================
            if not skip_modules.get('abricate', False):
                self.banner.display_module_header("ABRICATE ANALYSIS", "Comprehensive Resistance, Plasmid & Virulence Gene Screening")
                ab_success = self.run_abricate_analysis(fasta_files, output_path, threads,
                                                        min_identity=abricate_min_identity,
                                                        min_coverage=abricate_min_coverage)
                if ab_success:
                    self.banner.display_success("✅ ABRicate completed")
                else:
                    self.banner.display_error("❌ ABRicate failed")
                self.banner.display_random_quote()
                print()

            if not skip_modules.get('amr', False):
                self.banner.display_module_header("AMR ANALYSIS", "Antimicrobial Resistance Gene Detection with Point Mutations")
                amr_success = self.run_amr_analysis(fasta_files, output_path, threads,
                                                    min_identity=amr_min_identity,
                                                    min_coverage=amr_min_coverage,
                                                    skip_mutations=amr_skip_mutations,
                                                    force_update=amr_force_update)
                if amr_success:
                    self.banner.display_success("✅ AMR completed")
                else:
                    self.banner.display_error("❌ AMR failed")
                self.banner.display_random_quote()
                print()

            # ================================================================
            # Final summary
            # ================================================================
            if not skip_summary:
                self.banner.display_module_header("ULTIMATE REPORTER", "Gene‑centric Integrated Analysis")
                summary_success = self.run_summary_analysis(output_path)
                if summary_success:
                    self.banner.display_success("✅ Ultimate Reporter completed")
                else:
                    self.banner.display_warning("Ultimate Reporter had issues")
                self.banner.display_random_quote()
                print()

            # ================================================================
            # Final summary and citation
            # ================================================================
            analysis_time = datetime.now() - start_time
            self.banner.display_footer(analysis_time=str(analysis_time).split('.')[0],
                                       samples_processed=len(fasta_files))

            print(f"{Color.BOLD}{Color.BRIGHT_CYAN}📚 Please cite our EnteroScope paper:{Color.RESET}")
            print("   Beckley, B. EnteroScope: A Tailored Computational Workflow Enabling Rapid, User-Friendly Genotyping and Epidemiological Surveillance of the Enterobacter cloacae Complex. (In preparation)")
            print(f"{Color.BOLD}{Color.BRIGHT_CYAN}🔗 Support & Contributions:{Color.RESET}")
            print("   • Issues & feature requests: https://github.com/bbeckley-hub/enteroscope/issues")
            print("   • Email: brownbeckley94@gmail.com")
            print(f"{Color.YELLOW}⭐ Star us on GitHub if you find this tool useful! ⭐{Color.RESET}")

            # Final closing quote
            self.banner.display_random_quote()

        except KeyboardInterrupt:
            self.banner.display_error("Analysis interrupted by user")
        except Exception as e:
            self.banner.display_error(f"Critical error in analysis pipeline: {str(e)}")
            traceback.print_exc()


# =============================================================================
# Command-Line Entry Point
# =============================================================================
def main():
    """Parse command line arguments and run the EnteroScope orchestrator."""
    parser = argparse.ArgumentParser(
        description="EnteroScope: Advanced Enterobacter cloacae Complex Genomic Analysis Pipeline with FASTA QC & Ultimate Reporting",
        formatter_class=ColoredHelpFormatter,
        epilog=f"""
{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Examples:{ColoredHelpFormatter.RESET}
  # Basic analysis
  enteroscope -i genome.fna -o results/

  # Batch with glob pattern
  enteroscope -i "*.fna" -o batch_results --threads 8

  # Skip some modules
  enteroscope -i "*.fasta" -o analysis --threads 16 --skip-abricate --skip-amr

  # AMR with custom thresholds and no mutation reporting
  enteroscope -i "*.fna" -o results --amr-min-identity 0.95 --amr-min-coverage 0.9 --skip-amr-mutations

  # ABRicate with custom thresholds
  enteroscope -i "*.fna" -o results --abricate-min-identity 85 --abricate-min-coverage 85

  # Force update AMR database before analysis
  enteroscope -i "*.fna" -o results --amr-force-update

  # Update AMR database only (standalone)
  enteroscope --update-amr-db

  # Force update AMR database only
  enteroscope --force-update-amr-db

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Supported FASTA formats:{ColoredHelpFormatter.RESET} .fna, .fasta, .fa, .fn

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Analysis Modules:{ColoredHelpFormatter.RESET}
  • FASTA QC (Quality Control & Assembly Statistics)
  • MLST (Multi‑Locus Sequence Typing)
  • ABRicate (Comprehensive resistance, plasmid, virulence, BACMET screening)
  • AMR (AMRfinderPlus with point mutations – enabled by default)
  • Ultimate Reporter (Gene‑centric integrated HTML dashboard)

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Output:{ColoredHelpFormatter.RESET}
  Comprehensive results in organized subdirectories.
  A detailed log file (enteroscope_run.log) is written to the output directory.

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Citation:{ColoredHelpFormatter.RESET}
  Beckley B. EnteroScope: a species-optimized computational pipeline for rapid and accessible
  Enterobacter cloacae complex genotyping and surveillance. (In preparation)

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Support & Contributions:{ColoredHelpFormatter.RESET}
  • Issues & feature requests: https://github.com/bbeckley-hub/enteroscope/issues
  • Email: brownbeckley94@gmail.com

{ColoredHelpFormatter.YELLOW}⭐ Star us on GitHub if you find this tool useful! ⭐{ColoredHelpFormatter.RESET}
        """
    )
    # Required arguments
    parser.add_argument('-i', '--input', help='Input FASTA file(s) – can use glob patterns like "*.fna" or "*.fasta"')
    parser.add_argument('-o', '--output', help='Output directory for all results')
    parser.add_argument('-t', '--threads', type=int, default=2, help='Number of threads (default: 2)')

    # Skip options
    parser.add_argument('--skip-qc', action='store_true', help='Skip FASTA QC analysis')
    parser.add_argument('--skip-mlst', action='store_true', help='Skip MLST analysis')
    parser.add_argument('--skip-abricate', action='store_true', help='Skip ABRicate analysis')
    parser.add_argument('--skip-amr', action='store_true', help='Skip AMR analysis (AMRfinderPlus)')
    parser.add_argument('--skip-summary', action='store_true', help='Skip ultimate reporter generation')

    # AMR‑specific flags
    parser.add_argument('--amr-min-identity', type=float, help='Minimum identity (0..1) for AMR hits (default: auto)')
    parser.add_argument('--amr-min-coverage', type=float, help='Minimum coverage (0..1) for AMR hits (default: auto)')
    parser.add_argument('--skip-amr-mutations', action='store_true', help='Disable point mutation reporting in AMR (enabled by default)')
    parser.add_argument('--amr-force-update', action='store_true', help='Force update AMR database before analysis')

    # ABRicate‑specific flags
    parser.add_argument('--abricate-min-identity', type=float, help='Minimum identity percentage for ABRicate (default: 80)')
    parser.add_argument('--abricate-min-coverage', type=float, help='Minimum coverage percentage for ABRicate (default: 80)')

    # Database update standalone
    parser.add_argument('--update-amr-db', action='store_true', help='Update AMRfinderPlus database (incremental) and exit')
    parser.add_argument('--force-update-amr-db', action='store_true', help='Force complete AMR database update (overwrites old) and exit')

    # Utility flags
    parser.add_argument('--keep-temp', action='store_true', help='Do not delete temporary directories (for debugging)')
    parser.add_argument('--clean-output', action='store_true', help='Delete output directory before analysis (prevents mixing results from different runs)')
    parser.add_argument('--version', action='version', version=f'EnteroScope v{__version__}')

    args = parser.parse_args()

    # Standalone AMR database update
    if args.update_amr_db or args.force_update_amr_db:
        orch = EnteroScopeOrchestrator()
        if args.force_update_amr_db:
            orch.update_amr_database(force=True)
        else:
            orch.update_amr_database(force=False)
        sys.exit(0)

    # For analysis, input and output are required
    if not args.input or not args.output:
        parser.error("When not using --update-amr-db or --force-update-amr-db, both -i/--input and -o/--output are required.")

    skip_modules = {
        'qc': args.skip_qc,
        'mlst': args.skip_mlst,
        'abricate': args.skip_abricate,
        'amr': args.skip_amr,
    }

    orch = EnteroScopeOrchestrator()
    orch.keep_temp = args.keep_temp
    orch.run_complete_analysis(
        input_path=args.input,
        output_dir=args.output,
        threads=args.threads,
        skip_modules=skip_modules,
        skip_summary=args.skip_summary,
        amr_min_identity=args.amr_min_identity,
        amr_min_coverage=args.amr_min_coverage,
        amr_skip_mutations=args.skip_amr_mutations,
        amr_force_update=args.amr_force_update,
        abricate_min_identity=args.abricate_min_identity,
        abricate_min_coverage=args.abricate_min_coverage,
        clean_output=args.clean_output
    )


if __name__ == "__main__":
    main()