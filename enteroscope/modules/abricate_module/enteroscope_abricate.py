#!/usr/bin/env python3
"""
EnteroScope ABRicate Standalone Module
Comprehensive ABRicate analysis with HTML, TSV, and JSON reporting - MAXIMUM SPEED VERSION
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
GitHub: https://github.com/bbeckley-hub/enteroscope
Date: 2026
Version: 1.0.0 (teal/cyan theme)
"""

import subprocess
import sys
import os
import glob
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any, Set, Tuple
import argparse
import re
from datetime import datetime
import psutil
import math
import json
import random 
from collections import defaultdict, Counter

class EnteroAbricateExecutor:
    """ABRicate executor for Enterobacter cloacae complex - MAXIMUM SPEED (capped at 64 cores)"""
    
    def __init__(self, cpus: int = None):
        self.logger = self._setup_logging()
        self.available_ram = self._get_available_ram()
        self.cpus = self._calculate_optimal_cpus(cpus)

        # Databases suitable for Enterobacteriaceae
        self.required_databases = [
            'ncbi', 'card', 'resfinder', 'vfdb', 'argannot',
            'plasmidfinder', 'megares', 'ecoh', 'ecoli_vf', 'bacmet2'
        ]

        # ---------- Enterobacteriaceae / E. cloacae gene sets ----------
        # Critical carbapenemases
        self.critical_carbapenemases = {
            'blaKPC', 'blaNDM', 'blaVIM', 'blaIMP', 'blaOXA-48', 'blaOXA-181', 'blaOXA-232',
            'blaGES', 'blaIMI', 'blaSME', 'blaNMC', 'blaCcrA', 'blaCphA',
            'KPC', 'NDM', 'VIM', 'IMP', 'OXA-48', 'OXA-181', 'OXA-232',
            'GES', 'IMI', 'SME', 'NMC', 'CcrA', 'CphA'
        }

        # Critical ESBLs / AmpC
        self.critical_esbls = {
            'blaCTX-M', 'blaSHV', 'blaTEM', 'blaPER', 'blaVEB', 'blaGES', 'blaBEL',
            'CTX-M', 'SHV', 'TEM', 'PER', 'VEB', 'GES', 'BEL', 'BES',
            'blaCMY', 'blaDHA', 'blaFOX', 'blaMOX', 'blaACC', 'blaACT', 'blaMIR',
            'CMY', 'DHA', 'FOX', 'MOX', 'ACC', 'ACT', 'MIR'
        }

        # Colistin resistance
        self.critical_colistin = {
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10',
            'mcr1', 'mcr2', 'mcr3', 'mcr4', 'mcr5', 'mcr6', 'mcr7', 'mcr8', 'mcr9', 'mcr10',
            'pmrA', 'pmrB', 'pmrC', 'lpxA', 'lpxC', 'lpxD', 'arnA', 'arnB', 'arnC', 'arnD', 'eptA'
        }

        # Critical aminoglycoside resistance
        self.critical_aminoglycoside = {
            'armA', 'rmtA', 'rmtB', 'rmtC', 'rmtD', 'rmtE', 'rmtF', 'rmtG', 'rmtH',
            'APH(3\')', 'APH(3\')-VI', 'APH(6)', 'AAC(3)', 'AAC(6\')', 'ANT(2")', 'ANT(3")', 'ANT(4")',
            'aacC1', 'aacC2', 'aacC4', 'aacA4', 'aacA7', 'aadA1', 'aadA2', 'aadA5', 'aadA7',
            'strA', 'strB', 'aphA1', 'aphA2', 'aphA3', 'aphA6', 'aac3', 'aac6', 'aadA', 'aadB'
        }

        # High‑risk resistance (tetracycline, sulfonamide, trimethoprim, chloramphenicol, macrolide, quinolone, efflux pumps)
        self.high_risk_resistance = {
            'tetA', 'tetB', 'tetC', 'tetD', 'tetE', 'tetG', 'tetH', 'tetK', 'tetL', 'tetM', 'tetO', 'tetQ', 'tetS', 'tetX',
            'sul1', 'sul2', 'sul3', 'sul4',
            'dfrA1', 'dfrA5', 'dfrA7', 'dfrA8', 'dfrA12', 'dfrA14', 'dfrA17', 'dfrA19', 'dfrA20', 'dfrA21',
            'dfrB1', 'dfrB2', 'dfrB3', 'dfrB4', 'dfrB5', 'dfrB6', 'dfrB7',
            'catA1', 'catA2', 'catB2', 'catB3', 'catB8', 'catI', 'catII', 'catIII',
            'cmlA', 'cmlA1', 'cmlA5', 'cmlA6', 'cmlA7', 'floR',
            'ermA', 'ermB', 'ermC', 'ermF', 'ermG', 'ermX', 'ermY',
            'mphA', 'mphB', 'mphC', 'mphD', 'mphE', 'msrA', 'msrB', 'msrC', 'msrD',
            'qnrA1', 'qnrB1', 'qnrB2', 'qnrB4', 'qnrB6', 'qnrB10', 'qnrB19', 'qnrS1', 'qnrS2', 'qnrVC1', 'qnrVC4',
            'aac(6\')-Ib-cr', 'qepA', 'qepA1', 'qepA2', 'qepA3', 'qepA4',
            'fosA', 'fosB', 'fosC', 'fosX',
            'arr-2', 'arr-3', 'arr-4', 'arr-5', 'arr-6', 'arr-7',
            'adeA', 'adeB', 'adeC', 'adeG', 'adeH', 'adeI', 'adeJ', 'adeK', 'adeL', 'adeM', 'adeN',
            'abeM', 'abeS', 'acrA', 'acrB', 'tolC', 'mexA', 'mexB', 'mexC', 'mexD', 'mexE', 'mexF',
            'oprM', 'oprN', 'oprJ', 'mdtA', 'mdtB', 'mdtC', 'emrA', 'emrB', 'emrD', 'emrE', 'emrK', 'emrY',
            'mdx', 'acrAB', 'acrEF', 'mdfA', 'cmlA'
        }

        # Virulence factors (Enterobacteriaceae)
        self.virulence_genes = {
            'csuA', 'csuB', 'csuC', 'csuD', 'csuE', 'csuA/B', 'csuC/D/E',
            'bfmR', 'bfmS', 'ompA', 'bap', 'epsA', 'ptk',
            'basA', 'basB', 'basC', 'basD', 'basE', 'basF', 'basG', 'basH', 'basI', 'basJ',
            'barA', 'barB', 'bauA', 'bauB', 'bauC', 'bauD', 'bauE', 'bauF',
            'entA', 'entB', 'entC', 'entD', 'entE', 'entF', 'entS',
            'fhuA', 'fhuB', 'fhuC', 'fhuD', 'fhuE', 'fhuF',
            'fepA', 'fepB', 'fepC', 'fepD', 'fepE', 'fepG',
            'iroN', 'iucA', 'iucB', 'iucC', 'iucD', 'iutA',
            'abaI', 'abaR', 'luxI', 'luxR', 'qseC', 'qseB',
            'hlyA', 'hlyB', 'hlyC', 'hlyD', 'rtxA', 'rtxB', 'rtxC', 'rtxD', 'rtxE',
            'tssA', 'tssB', 'tssC', 'tssD', 'tssE', 'tssF', 'tssG', 'tssH', 'tssI', 'tssJ', 'tssK', 'tssL', 'tssM',
            'vgrG', 'vgrG1', 'vgrG2', 'vgrG3', 'vgrG4', 'vgrG5',
            'cvaA', 'cvaB', 'cvaC', 'colicin', 'microcin',
            'plc1', 'plc2', 'plcD', 'lipA', 'lipB', 'lipC', 'lipH',
            'aprA', 'aprB', 'aprC', 'aprD', 'aprE', 'aprF',
            'cpsA', 'cpsB', 'cpsC', 'cpsD', 'cpsE', 'cpsF', 'cpsG', 'cpsH', 'cpsI', 'cpsJ',
            'kpsM', 'kpsS', 'kpsT', 'kpsU', 'kpsV', 'kpsW',
            'wza', 'wzb', 'wzc', 'wzd', 'wze', 'wzf', 'wzg',
            'fimA', 'fimB', 'fimC', 'fimD', 'fimE', 'fimF', 'fimG', 'fimH',
            'papA', 'papB', 'papC', 'papD', 'papE', 'papF', 'papG', 'papH',
            'sfmA', 'sfmB', 'sfmC', 'sfmD', 'sfmE', 'sfmF', 'sfmG', 'sfmH',
            'type1', 'P-pili', 'curli', 'csgA', 'csgB', 'csgD',
            'pilA', 'pilB', 'pilC', 'pilD', 'pilE', 'pilF', 'pilG', 'pilH', 'pilI', 'pilJ', 'pilK', 'pilL', 'pilM', 'pilN', 'pilO', 'pilP', 'pilQ'
        }

        # Combined risk categories
        self.critical_resistance_genes = (
            self.critical_carbapenemases.union(self.critical_esbls)
            .union(self.critical_colistin).union(self.critical_aminoglycoside)
        )
        self.critical_virulence_genes = {
            'hlyA', 'hlyB', 'hlyC', 'hlyD', 'rtxA', 'rtxB', 'rtxC', 'rtxD', 'rtxE',
            'cvaA', 'cvaB', 'cvaC', 'colicin', 'microcin'
        }

        self.metadata = {
            "tool_name": "EnteroScope ABRicate",
            "version": "1.0.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub/enteroscope",
            "affiliation": "University of Ghana Medical School - Department of Medical Biochemistry",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }

        self.science_quotes = [
            {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein"},
            {"text": "Science is not only a disciple of reason but also one of romance and passion.", "author": "Stephen Hawking"},
            {"text": "Somewhere, something incredible is waiting to be known.", "author": "Carl Sagan"},
            {"text": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson"},
            {"text": "In science, there are no shortcuts to truth.", "author": "Karl Popper"},
            {"text": "Science knows no country, because knowledge belongs to humanity.", "author": "Louis Pasteur"},
            {"text": "The science of today is the technology of tomorrow.", "author": "Edward Teller"},
            {"text": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie"},
            {"text": "Research is what I'm doing when I don't know what I'm doing.", "author": "Wernher von Braun"},
            {"text": "The universe is not required to be in perfect harmony with human ambition.", "author": "Carl Sagan"},
            {"text": "EnteroScope represents the convergence of genomic surveillance and clinical diagnostics for Enterobacter cloacae complex.", "author": "Brown Beckley"},
            {"text": "In the battle against antimicrobial resistance, EnteroScope reveals the genetic blueprints of resistant Enterobacter strains.", "author": "Brown Beckley"},
            {"text": "EnteroScope bridges the gap between sequencing data and public health action for ECC.", "author": "Brown Beckley"},
            {"text": "Through EnteroScope, we turn complexity into clear reports, empowering clinicians and researchers.", "author": "Brown Beckley"},
            {"text": "EnteroScope makes advanced genomic typing for Enterobacter cloacae complex accessible to all.", "author": "Brown Beckley"}
        ]

        # ASCII art for EnteroScope 
        self.ascii_art = r"""
███████╗███╗   ██╗████████╗███████╗██████╗  ██████╗ ███████╗ ██████╗  ██████╗ ██████╗ ███████╗
██╔════╝████╗  ██║╚══██╔══╝██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔════╝ ██╔═══██╗██╔══██╗██╔════╝
█████╗  ██╔██╗ ██║   ██║   █████╗  ██████╔╝██║   ██║███████╗██║      ██║   ██║██████╔╝█████╗  
██╔══╝  ██║╚██╗██║   ██║   ██╔══╝  ██╔══██╗██║   ██║╚════██║██║      ██║   ██║██╔═══╝ ██╔══╝  
███████╗██║ ╚████║   ██║   ███████╗██║  ██║╚██████╔╝███████║╚██████╗ ╚██████╔╝██║     ███████╗
╚══════╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝ ╚═════╝  ╚═════╝ ╚═╝     ╚══════╝
"""

    def get_random_quote(self):
        return random.choice(self.science_quotes)

    def _setup_logging(self):
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        return logging.getLogger(__name__)

    def _get_available_ram(self) -> int:
        try:
            return psutil.virtual_memory().available / (1024 ** 3)
        except Exception:
            return 8

    def _calculate_optimal_cpus(self, user_cpus: int = None) -> int:
        if user_cpus is not None:
            user_cpus = min(user_cpus, 64)
            self._log_resource_info(user_cpus)
            return user_cpus
        try:
            total = psutil.cpu_count(logical=False) or os.cpu_count() or 2
            total = min(total, 64)
            if total <= 4:
                opt = total
            elif total <= 8:
                opt = total - 1
            elif total <= 16:
                opt = max(8, total - 2)
            elif total <= 32:
                opt = max(16, total - 3)
            else:
                opt = min(64, int(total * 0.95))
            opt = max(1, min(opt, total, 64))
            self._log_resource_info(opt, total)
            return opt
        except Exception:
            return min(os.cpu_count() or 4, 64)

    def _log_resource_info(self, cpus: int, total_cores: int = None):
        self.logger.info(f"Available RAM: {self.available_ram:.1f} GB")
        if total_cores:
            self.logger.info(f"System CPU cores: {total_cores}")
            self.logger.info(f"Using CPU cores: {cpus} ({cpus/total_cores*100:.1f}%)")
        else:
            self.logger.info(f"Using user‑specified CPU cores: {cpus}")
        if cpus >= 48:
            self.logger.info("💡 Performance: EXTREME SPEED MODE 🚀 (>48 cores)")
        elif cpus >= 32:
            self.logger.info("💡 Performance: ULTRA MAXIMUM SPEED MODE 🚀")
        elif cpus >= 16:
            self.logger.info("💡 Performance: MAXIMUM SPEED MODE 🚀")
        elif cpus >= 8:
            self.logger.info("💡 Performance: High‑speed mode")
        else:
            self.logger.info("💡 Performance: Standard multi‑core")

    def check_abricate_installed(self) -> bool:
        try:
            subprocess.run(['abricate', '--version'], capture_output=True, text=True, check=True)
            self.logger.info("ABRicate is installed.")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.logger.error("ABRicate not found. Install with: conda install -c bioconda abricate")
            return False

    def setup_abricate_databases(self):
        self.logger.info("Setting up ABRicate databases for Enterobacter analysis...")
        available, missing = [], []
        try:
            out = subprocess.run(['abricate', '--list'], capture_output=True, text=True, check=True).stdout
            for db in self.required_databases:
                if db in out:
                    self.logger.info(f"✓ Database available: {db}")
                    available.append(db)
                else:
                    missing.append(db)
            for db in missing:
                self.logger.info(f"Setting up {db}...")
                subprocess.run(['abricate', '--setupdb', '--db', db], capture_output=True, text=True, check=True)
                self.logger.info(f"✓ Database setup completed: {db}")
                available.append(db)
            self.required_databases = available
            self.logger.info(f"Using databases: {', '.join(self.required_databases)}")
        except Exception as e:
            self.logger.error(f"Database setup error: {e}")

    def run_abricate_single_db(self, genome_file: str, database: str, out_dir: str) -> Dict:
        genome_name = Path(genome_file).stem
        out_file = os.path.join(out_dir, f"abricate_{database}.txt")
        cmd = ['abricate', genome_file, '--db', database, '--minid', '80', '--mincov', '80']
        self.logger.info(f"Running ABRicate: {genome_name} --db {database}")
        try:
            with open(out_file, 'w') as f:
                subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True, check=True)
            hits = self._parse_abricate_output(out_file)
            self._create_database_html_report(genome_name, database, hits, out_dir)
            return {'database': database, 'genome': genome_name, 'output_file': out_file,
                    'hits': hits, 'hit_count': len(hits), 'status': 'success'}
        except subprocess.CalledProcessError as e:
            self.logger.error(f"ABRicate failed for {database} on {genome_name}: {e.stderr}")
            return {'database': database, 'genome': genome_name, 'output_file': out_file,
                    'hits': [], 'hit_count': 0, 'status': 'failed'}

    def _parse_abricate_output(self, abricate_file: str) -> List[Dict]:
        hits = []
        try:
            with open(abricate_file, 'r') as f:
                lines = f.readlines()
            if not lines:
                return hits
            headers, data_lines = [], []
            for line in lines:
                if line.startswith('#FILE') and not headers:
                    headers = line.strip().replace('#', '').split('\t')
                elif line.strip() and not line.startswith('#'):
                    data_lines.append(line.strip())
            if not headers:
                return hits
            expected = len(headers)
            for line in data_lines:
                parts = line.split('\t')
                if len(parts) > expected:
                    parts = parts[:expected-1] + ['\t'.join(parts[expected-1:])]
                elif len(parts) < expected:
                    parts += [''] * (expected - len(parts))
                if len(parts) == expected:
                    hit = {headers[i]: parts[i] for i in range(expected)}
                    processed = {
                        'file': hit.get('FILE', ''),
                        'sequence': hit.get('SEQUENCE', ''),
                        'start': hit.get('START', ''),
                        'end': hit.get('END', ''),
                        'strand': hit.get('STRAND', ''),
                        'gene': hit.get('GENE', ''),
                        'coverage': hit.get('COVERAGE', ''),
                        'coverage_map': hit.get('COVERAGE_MAP', ''),
                        'gaps': hit.get('GAPS', ''),
                        'coverage_percent': hit.get('%COVERAGE', ''),
                        'identity_percent': hit.get('%IDENTITY', ''),
                        'database': hit.get('DATABASE', ''),
                        'accession': hit.get('ACCESSION', ''),
                        'product': hit.get('PRODUCT', ''),
                        'resistance': hit.get('RESISTANCE', '')
                    }
                    hits.append(processed)
        except Exception as e:
            self.logger.error(f"Error parsing {abricate_file}: {e}")
        return hits

    # ---------- Individual database HTML report (teal/cyan) ----------
    def _create_database_html_report(self, genome_name: str, database: str, hits: List[Dict], out_dir: str):
        random_quote = self.get_random_quote()
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EnteroScope ABRicate - {database.upper()} Database Report</title>
    <style>
        * {{ margin:0; padding:0; box-sizing:border-box; }}
        body {{
            background: linear-gradient(135deg, #0f766e 0%, #0d9488 50%, #06b6d4 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{
            background: rgba(0,0,0,0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            border: 2px solid #06b6d4;
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            white-space: pre;
            color: #06b6d4;
            text-shadow: 0 0 10px #06b6d4;
            overflow-x: auto;
        }}
        .quote-container {{
            background: rgba(0,0,0,0.4);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            border-left: 4px solid #06b6d4;
            transition: opacity 0.5s ease-in-out;
        }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; color: #fff; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
        .report-section {{
            background: rgba(255,255,255,0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 12px;
            margin-bottom: 20px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.2);
        }}
        .report-section h2 {{
            color: #0f766e;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
        }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-top: 15px;
        }}
        .metric-card {{
            background: linear-gradient(135deg, #06b6d4 0%, #0891b2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.15);
        }}
        .metric-label {{ font-size: 14px; opacity: 0.9; margin-bottom: 5px; }}
        .metric-value {{ font-size: 24px; font-weight: bold; }}
        .summary-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            font-size: 14px;
        }}
        .summary-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #2563eb 100%);
            color: white;
            padding: 12px;
            text-align: left;
        }}
        .summary-table td {{
            padding: 12px;
            border-bottom: 1px solid #e5e7eb;
        }}
        .summary-table tr:nth-child(even) {{ background-color: #f8fafc; }}
        .summary-table tr:hover {{ background-color: #e0f2fe; }}
        .table-responsive {{ overflow-x: auto; margin: 20px 0; }}
        .present {{ background-color: #d1fae5; }}
        .high-risk {{ background-color: #fef3c7; }}
        .critical {{ background-color: #fee2e2; font-weight: bold; }}
        .product-cell {{
            white-space: normal !important;
            word-wrap: break-word;
            max-width: 400px;
            min-width: 200px;
        }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0,0,0,0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        .footer a {{ color: #fbbf24; text-decoration: none; }}
        .footer a:hover {{ text-decoration: underline; }}
        .timestamp {{ color: #fbbf24; font-weight: bold; }}
        .authorship {{
            margin-top: 15px;
            padding: 15px;
            background: rgba(255,255,255,0.1);
            border-radius: 8px;
            font-size: 12px;
        }}
        @media (max-width: 768px) {{ .ascii-art {{ font-size: 6px; }} }}
    </style>
</head>
<body>
<div class="container">
    <div class="header">
        <div class="ascii-container"><div class="ascii-art">{self.ascii_art}</div></div>
        <div class="quote-container" id="quoteContainer">
            <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
            <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
        </div>
    </div>
    <div class="report-section">
        <h2>📊 Database Information</h2>
        <div class="metrics-grid">
            <div class="metric-card"><div class="metric-label">Database</div><div class="metric-value">{database.upper()}</div></div>
            <div class="metric-card"><div class="metric-label">Genome</div><div class="metric-value">{genome_name}</div></div>
            <div class="metric-card"><div class="metric-label">Analysis Date</div><div class="metric-value">{datetime.now().strftime('%Y-%m-%d %H:%M')}</div></div>
            <div class="metric-card"><div class="metric-label">Total Hits</div><div class="metric-value">{len(hits)}</div></div>
        </div>
    </div>
"""
        if hits:
            html += """
    <div class="report-section">
        <h2>🔍 Genes Detected</h2>
        <div class="table-responsive">
            <table class="summary-table">
                <thead><tr><th>Gene</th><th>Product</th><th>Coverage</th><th>Identity</th><th>Accession</th></tr></thead>
                <tbody>
"""
            for hit in hits:
                row_class = "critical" if any(c in hit['gene'] for c in self.critical_resistance_genes) else ("high-risk" if any(v in hit['gene'] for v in self.virulence_genes) else "present")
                html += f'<tr class="{row_class}"><td><strong>{hit["gene"]}</strong></td><td class="product-cell">{hit["product"]}</td><td>{hit["coverage_percent"]}%</td><td>{hit["identity_percent"]}%</td><td>{hit["accession"]}</td></tr>'
            html += """
                </tbody>
            </table>
        </div>
    </div>
"""
        else:
            html += '<div class="report-section"><p>No significant hits found.</p></div>'

        html += f"""
    <div class="footer">
        <p><strong>EnteroScope</strong> - ABRicate Analysis Module</p>
        <p class="timestamp">Generated: {current_time}</p>
        <div class="authorship">
            <p><strong>Author:</strong> Brown Beckley | <strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub/enteroscope" target="_blank">https://github.com/bbeckley-hub/enteroscope</a></p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
        </div>
    </div>
</div>
<script>
    const quotes = {json.dumps(self.science_quotes)};
    const container = document.getElementById('quoteContainer');
    const textDiv = document.getElementById('quoteText');
    const authorDiv = document.getElementById('quoteAuthor');
    let idx = 0;
    function rotate() {{
        container.style.opacity = '0';
        setTimeout(() => {{
            const q = quotes[idx];
            textDiv.innerHTML = '"' + q.text + '"';
            authorDiv.innerHTML = '— ' + q.author;
            container.style.opacity = '1';
            idx = (idx + 1) % quotes.length;
        }}, 500);
    }}
    setInterval(rotate, 10000);
    document.addEventListener('DOMContentLoaded', rotate);
</script>
</body>
</html>"""
        out_file = os.path.join(out_dir, f"abricate_{database}_report.html")
        with open(out_file, 'w') as f:
            f.write(html)
        self.logger.info(f"Individual database report: {out_file}")

    # ---------- Resistance analysis for Enterobacter ----------
    def analyze_enterobacter_resistance(self, all_hits: List[Dict]) -> Dict[str, Any]:
        analysis = {
            'carbapenemase_status': 'negative',
            'esbl_status': 'negative',
            'colistin_resistance': 'negative',
            'critical_carbapenemase_genes': [],
            'critical_esbl_genes': [],
            'critical_colistin_genes': [],
            'critical_aminoglycoside_genes': [],
            'high_risk_resistance_genes': [],
            'critical_virulence_genes': [],
            'high_risk_virulence_genes': [],
            'moderate_risk_genes': [],
            'resistance_classes': {},
            'total_critical_resistance': 0,
            'total_high_risk_resistance': 0,
            'total_critical_virulence': 0,
            'total_high_risk_virulence': 0,
            'total_hits': len(all_hits)
        }
        for hit in all_hits:
            gene = hit['gene']
            gl = gene.lower()
            # Carbapenemase
            if any(x in gl for x in ['kpc', 'ndm', 'vim', 'imp', 'oxa-48', 'ges', 'imi', 'sme', 'nmc']):
                if any(c in gl for c in [x.lower() for x in self.critical_carbapenemases]):
                    analysis['carbapenemase_status'] = 'positive'
                    analysis['critical_carbapenemase_genes'].append(self._gene_dict(hit, 'CARBAPENEMASE'))
            # ESBL
            elif any(x in gl for x in ['ctx-m', 'shv', 'tem', 'per', 'veb', 'ges', 'cmy', 'dha', 'fox']):
                if any(c in gl for c in [x.lower() for x in self.critical_esbls]):
                    analysis['esbl_status'] = 'positive'
                    analysis['critical_esbl_genes'].append(self._gene_dict(hit, 'ESBL'))
            # Colistin
            elif 'mcr' in gl:
                if any(c in gl for c in [x.lower() for x in self.critical_colistin]):
                    analysis['colistin_resistance'] = 'positive'
                    analysis['critical_colistin_genes'].append(self._gene_dict(hit, 'COLISTIN-RES'))
            # Aminoglycoside
            elif any(x in gl for x in ['arm', 'rmt', 'aph', 'aac', 'aad', 'str', 'ant']):
                if any(c in gl for c in [x.lower() for x in self.critical_aminoglycoside]):
                    analysis['critical_aminoglycoside_genes'].append(self._gene_dict(hit, 'CRITICAL_AMINOGLYCOSIDE'))
            # High‑risk resistance
            elif any(hr in gl for hr in [x.lower() for x in self.high_risk_resistance]):
                analysis['high_risk_resistance_genes'].append(self._gene_dict(hit, 'HIGH-RISK'))
            # Virulence
            elif any(vf in gl for vf in [x.lower() for x in self.virulence_genes]):
                if any(cv in gl for cv in ['hly', 'rtx', 'colicin', 'microcin']):
                    analysis['critical_virulence_genes'].append(self._gene_dict(hit, 'CRITICAL-VIRULENCE'))
                else:
                    analysis['high_risk_virulence_genes'].append(self._gene_dict(hit, 'HIGH-VIRULENCE'))
            # Resistance class
            rclass = self._classify_resistance(hit['product'], gene)
            if rclass:
                analysis['resistance_classes'].setdefault(rclass, [])
                if not any(g['gene'] == gene for g in analysis['resistance_classes'][rclass]):
                    analysis['resistance_classes'][rclass].append({'gene': gene, 'product': hit['product']})
        analysis['total_critical_resistance'] = sum(len(analysis[k]) for k in ['critical_carbapenemase_genes', 'critical_esbl_genes', 'critical_colistin_genes', 'critical_aminoglycoside_genes'])
        analysis['total_high_risk_resistance'] = len(analysis['high_risk_resistance_genes'])
        analysis['total_critical_virulence'] = len(analysis['critical_virulence_genes'])
        analysis['total_high_risk_virulence'] = len(analysis['high_risk_virulence_genes'])
        return analysis

    def _gene_dict(self, hit: Dict, level: str) -> Dict:
        return {
            'gene': hit['gene'],
            'product': hit['product'],
            'database': hit['database'],
            'coverage': hit['coverage_percent'],
            'identity': hit['identity_percent'],
            'risk_level': level
        }

    def _classify_resistance(self, product: str, gene: str) -> str:
        p = product.lower()
        g = gene.lower()
        if any(t in p or t in g for t in ['carbapenem', 'kpc', 'ndm', 'vim', 'imp', 'oxa-48']):
            return 'Carbapenem resistance'
        if any(t in p or t in g for t in ['beta-lactam', 'esbl', 'ctx', 'shv', 'tem', 'per', 'veb', 'ges', 'cmy', 'dha']):
            return 'Beta-lactam resistance'
        if any(t in p or t in g for t in ['aminoglycoside', 'aac', 'aad', 'aph', 'str', 'arm', 'rmt']):
            return 'Aminoglycoside resistance'
        if any(t in p or t in g for t in ['tetracycline', 'tet']):
            return 'Tetracycline resistance'
        if any(t in p or t in g for t in ['sulfonamide', 'sul']):
            return 'Sulfonamide resistance'
        if any(t in p or t in g for t in ['trimethoprim', 'dfr']):
            return 'Trimethoprim resistance'
        if any(t in p or t in g for t in ['chloramphenicol', 'cat', 'cml', 'flor']):
            return 'Chloramphenicol resistance'
        if any(t in p or t in g for t in ['macrolide', 'erm', 'mph', 'msr']):
            return 'Macrolide resistance'
        if any(t in p or t in g for t in ['quinolone', 'qnr', 'qep']):
            return 'Quinolone resistance'
        if any(t in p or t in g for t in ['colistin', 'polymyxin', 'mcr']):
            return 'Polymyxin resistance'
        if any(t in p or t in g for t in ['fosfomycin', 'fos']):
            return 'Fosfomycin resistance'
        if any(t in p or t in g for t in ['rifampicin', 'arr']):
            return 'Rifampicin resistance'
        if any(t in p or t in g for t in ['efflux', 'ade', 'abe', 'mex', 'acr', 'emr']):
            return 'Efflux pumps'
        if any(t in p or t in g for t in ['virulence', 'toxin', 'hemolysin', 'siderophore']):
            return 'Virulence factors'
        return 'Other resistance'

    # ---------- Comprehensive HTML report (Enterobacter) ----------
    def create_comprehensive_html_report(self, genome_name: str, results: Dict, out_dir: str):
        all_hits = [h for db_res in results.values() for h in db_res['hits']]
        analysis = self.analyze_enterobacter_resistance(all_hits)
        random_quote = self.get_random_quote()
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EnteroScope ABRicate - Comprehensive Report</title>
    <style>
        * {{ margin:0; padding:0; box-sizing:border-box; }}
        body {{
            background: linear-gradient(135deg, #0f766e 0%, #0d9488 50%, #06b6d4 100%);
            font-family: 'Segoe UI', sans-serif;
            color: #ffffff;
            padding: 20px;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .ascii-container {{ background: rgba(0,0,0,0.7); padding:20px; border-radius:15px; margin-bottom:20px; border:2px solid #06b6d4; }}
        .ascii-art {{ font-family:'Courier New',monospace; font-size:10px; white-space:pre; color:#06b6d4; text-shadow:0 0 10px #06b6d4; overflow-x:auto; }}
        .quote-container {{ background: rgba(0,0,0,0.4); backdrop-filter:blur(10px); padding:20px; border-radius:10px; margin-bottom:30px; text-align:center; border-left:4px solid #06b6d4; }}
        .quote-text {{ font-size:18px; font-style:italic; margin-bottom:10px; }}
        .quote-author {{ font-size:14px; color:#fbbf24; font-weight:bold; }}
        .report-section {{ background: rgba(255,255,255,0.95); color:#1f2937; padding:25px; border-radius:12px; margin-bottom:20px; }}
        .report-section h2 {{ color:#0f766e; border-bottom:3px solid #3b82f6; padding-bottom:10px; margin-bottom:20px; }}
        .metrics-grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(200px,1fr)); gap:20px; margin:20px 0; }}
        .metric-card {{ background:linear-gradient(135deg,#06b6d4 0%,#0891b2 100%); color:white; padding:20px; border-radius:8px; text-align:center; }}
        .summary-table {{ width:100%; border-collapse:collapse; margin-top:20px; font-size:14px; }}
        .summary-table th {{ background:linear-gradient(135deg,#3b82f6 0%,#2563eb 100%); color:white; padding:12px; text-align:left; }}
        .summary-table td {{ padding:12px; border-bottom:1px solid #e5e7eb; }}
        .summary-table tr:nth-child(even) {{ background-color:#f8fafc; }}
        .summary-table tr:hover {{ background-color:#e0f2fe; }}
        .table-responsive {{ overflow-x:auto; margin:20px 0; }}
        .risk-badge {{ display:inline-block; background:#dc3545; color:white; padding:5px 10px; border-radius:15px; margin:2px; font-size:0.9em; }}
        .warning-badge {{ background:#f59e0b; color:black; }}
        .footer {{ text-align:center; margin-top:30px; padding:20px; background:rgba(0,0,0,0.3); border-radius:10px; }}
        .footer a {{ color:#fbbf24; text-decoration:none; }}
        .footer a:hover {{ text-decoration:underline; }}
        .timestamp {{ color:#fbbf24; font-weight:bold; }}
        .authorship {{ margin-top:15px; padding:15px; background:rgba(255,255,255,0.1); border-radius:8px; font-size:12px; }}
        .present {{ background-color:#d1fae5; }}
        .high-risk {{ background-color:#fef3c7; }}
        .critical {{ background-color:#fee2e2; font-weight:bold; }}
    </style>
</head>
<body>
<div class="container">
    <div class="ascii-container"><div class="ascii-art">{self.ascii_art}</div></div>
    <div class="quote-container"><div id="science-quote">"{random_quote['text']}"</div><div id="science-author">— {random_quote['author']}</div></div>
    <div class="report-section">
        <h2>📊 AMR Summary - {genome_name}</h2>
        <div class="metrics-grid">
            <div class="metric-card"><div class="metric-label">Total Genes</div><div class="metric-value">{analysis['total_hits']}</div></div>
            <div class="metric-card"><div class="metric-label">Critical Resistance</div><div class="metric-value">{analysis['total_critical_resistance']}</div></div>
            <div class="metric-card"><div class="metric-label">Virulence Factors</div><div class="metric-value">{analysis['total_high_risk_virulence'] + analysis['total_critical_virulence']}</div></div>
        </div>
    </div>
"""
        # Critical alerts
        alerts = []
        if analysis['carbapenemase_status'] == 'positive': alerts.append("🔴 CARBAPENEMASE DETECTED")
        if analysis['esbl_status'] == 'positive': alerts.append("🟡 ESBL DETECTED")
        if analysis['colistin_resistance'] == 'positive': alerts.append("🔴 COLISTIN RESISTANCE DETECTED")
        if alerts:
            html += '<div class="report-section" style="border-left:4px solid #dc3545;"><h2 style="color:#dc3545;">⚠️ CRITICAL RESISTANCE ALERTS</h2><div>'
            for a in alerts: html += f'<span class="risk-badge">{a}</span>'
            html += '</div></div>'

        # Resistance classes
        if analysis['resistance_classes']:
            html += '<div class="report-section"><h2>🧪 Resistance Classes Detected</h2>'
            for cls, genes in analysis['resistance_classes'].items():
                glist = ', '.join(g['gene'] for g in genes)
                html += f'<div><strong>{cls}</strong> ({len(genes)} genes)<br><span style="color:#666;">{glist}</span></div><br>'
            html += '</div>'

        # Critical tables
        for title, key in [("🔴 Critical Carbapenemase Genes", 'critical_carbapenemase_genes'),
                           ("🟡 Critical ESBL Genes", 'critical_esbl_genes'),
                           ("🔴 Colistin Resistance Genes", 'critical_colistin_genes'),
                           ("🟡 High‑Risk Resistance Genes", 'high_risk_resistance_genes')]:
            if analysis[key]:
                html += f'<div class="report-section"><h2>{title}</h2><div class="table-responsive"><table class="summary-table"><thead><tr><th>Gene</th><th>Product</th><th>Database</th><th>Coverage</th><th>Identity</th><th>Risk</th></tr></thead><tbody>'
                for g in analysis[key]:
                    badge_class = "risk-badge" if "CARBAPENEMASE" in g['risk_level'] or "COLISTIN" in g['risk_level'] else "warning-badge"
                    html += f'<tr class="critical"><td><strong>{g["gene"]}</strong></td><td>{g["product"][:100]}</td><td>{g["database"]}</td><td>{g["coverage"]}%</td><td>{g["identity"]}%</td><td><span class="{badge_class}">{g["risk_level"]}</span></td></tr>'
                html += '</tbody></table></div></div>'

        # Database summary
        html += '<div class="report-section"><h2>🗃️ Database Results</h2><div class="table-responsive"><table class="summary-table"><thead><tr><th>Database</th><th>Hits</th><th>Status</th></tr></thead><tbody>'
        for db, res in results.items():
            status = "✅ success" if res['status'] == 'success' else "❌ failed"
            html += f'<tr class="present"><td>{db}</td><td>{res["hit_count"]}</td><td>{status}</td></tr>'
        html += '</tbody></table></div></div>'

        html += f"""
    <div class="footer">
        <p><strong>EnteroScope</strong> - ABRicate Analysis Module</p>
        <p class="timestamp">Generated: {current_time}</p>
        <div class="authorship">
            <p><strong>Author:</strong> Brown Beckley | <strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub/enteroscope" target="_blank">https://github.com/bbeckley-hub/enteroscope</a></p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
        </div>
    </div>
</div>
<script>
    const quotes = {json.dumps(self.science_quotes)};
    let idx = 0;
    const quoteDiv = document.getElementById('science-quote');
    const authorSpan = document.getElementById('science-author');
    function rotate() {{
        const q = quotes[idx];
        quoteDiv.innerHTML = '"' + q.text + '"';
        authorSpan.innerHTML = '— ' + q.author;
        idx = (idx + 1) % quotes.length;
    }}
    setInterval(rotate, 10000);
    document.addEventListener('DOMContentLoaded', rotate);
</script>
</body>
</html>"""
        out_file = os.path.join(out_dir, f"{genome_name}_comprehensive_abricate_report.html")
        with open(out_file, 'w') as f:
            f.write(html)
        self.logger.info(f"Comprehensive HTML report: {out_file}")

    # ---------- Batch summary methods (TSV, JSON, HTML) ----------
    def create_database_summaries(self, all_results: Dict, out_base: str):
        self.logger.info("Creating database summary files...")
        db_hits = defaultdict(list)
        for gname, gres in all_results.items():
            for db, db_res in gres['results'].items():
                for hit in db_res['hits']:
                    h = hit.copy()
                    h['genome'] = gname
                    db_hits[db].append(h)
        for db, hits in db_hits.items():
            if hits:
                # TSV
                tsv_file = os.path.join(out_base, f"enteroscope_{db}_abricate_summary.tsv")
                headers = list(hits[0].keys())
                with open(tsv_file, 'w') as f:
                    f.write('\t'.join(headers) + '\n')
                    for h in hits:
                        f.write('\t'.join(str(h.get(k, '')) for k in headers) + '\n')
                self.logger.info(f"TSV summary: {tsv_file} ({len(hits)} hits)")
                # JSON
                self._create_db_json_summary(db, hits, out_base)
                # HTML batch summary (StaphScope style)
                self._create_db_html_summary(db, hits, out_base)

    def _create_db_json_summary(self, db: str, hits: List[Dict], out_base: str):
        genomes = list(set(h['genome'] for h in hits))
        gene_freq = Counter(h['gene'] for h in hits)
        summary = {
            "database": db,
            "total_hits": len(hits),
            "unique_genomes": len(genomes),
            "unique_genes": len(gene_freq),
            "most_frequent_genes": gene_freq.most_common(20),
            "per_genome": {gen: len([h for h in hits if h['genome'] == gen]) for gen in genomes}
        }
        out_file = os.path.join(out_base, f"enteroscope_{db}_summary.json")
        with open(out_file, 'w') as f:
            json.dump(summary, f, indent=2)
        self.logger.info(f"JSON summary: {out_file}")

    # ----- BATCH HTML SUMMARY (StaphScope style, adapted to teal/cyan) -----
    def _create_db_html_summary(self, database: str, hits: List[Dict], out_base: str):
        """Create HTML summary report for a specific database across all genomes - StaphScope style."""
        random_quote = self.get_random_quote()
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        unique_genomes = sorted(set(hit['genome'] for hit in hits))
        genes_per_genome = {}
        for hit in hits:
            genome = hit['genome']
            if genome not in genes_per_genome:
                genes_per_genome[genome] = set()
            genes_per_genome[genome].add(hit['gene'])

        # Gene frequency
        gene_freq = Counter(hit['gene'] for hit in hits)

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EnteroScope - Batch ABRicate Summary ({database.upper()})</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            background: linear-gradient(135deg, #0f766e 0%, #0d9488 50%, #06b6d4 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{
            max-width: 1800px;
            margin: 0 auto;
        }}
        .header {{
            text-align: center;
            margin-bottom: 30px;
        }}
        .ascii-container {{
            background: rgba(0, 0, 0, 0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            border: 2px solid #06b6d4;
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            white-space: pre;
            color: #06b6d4;
            text-shadow: 0 0 10px #06b6d4;
            overflow-x: auto;
        }}
        .quote-container {{
            background: rgba(0, 0, 0, 0.4);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            border-left: 4px solid #06b6d4;
            transition: opacity 0.5s ease-in-out;
        }}
        .quote-text {{
            font-size: 18px;
            font-style: italic;
            margin-bottom: 10px;
            color: #ffffff;
        }}
        .quote-author {{
            font-size: 14px;
            color: #fbbf24;
            font-weight: bold;
        }}
        .report-section {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 12px;
            margin-bottom: 20px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.2);
        }}
        .report-section h2 {{
            color: #0f766e;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
        }}
        .summary-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            font-size: 14px;
        }}
        .summary-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #2563eb 100%);
            color: white;
            padding: 12px;
            text-align: left;
            position: sticky;
            top: 0;
        }}
        .summary-table td {{
            padding: 10px;
            border-bottom: 1px solid #e5e7eb;
        }}
        .summary-table tr:nth-child(even) {{
            background-color: #f8fafc;
        }}
        .summary-table tr:hover {{
            background-color: #e0f2fe;
        }}
        .table-responsive {{
            overflow-x: auto;
            margin: 20px 0;
        }}
        .spa-type-cell {{
            font-weight: bold;
            color: #0f766e;
        }}
        .repeat-cell {{
            font-family: 'Courier New', monospace;
            background-color: #f0fdfa;
            color: #0f766e;
            font-weight: bold;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-bottom: 20px;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #06b6d4 0%, #0891b2 100%);
            color: white;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }}
        .stat-value {{
            font-size: 24px;
            font-weight: bold;
            margin-bottom: 5px;
        }}
        .stat-label {{
            font-size: 12px;
            opacity: 0.9;
        }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        .footer a {{
            color: #fbbf24;
            text-decoration: none;
        }}
        .footer a:hover {{
            text-decoration: underline;
        }}
        .timestamp {{
            color: #fbbf24;
            font-weight: bold;
        }}
        .authorship {{
            margin-top: 15px;
            padding: 15px;
            background: rgba(255, 255, 255, 0.1);
            border-radius: 8px;
            font-size: 12px;
        }}
        @media (max-width: 768px) {{
            .ascii-art {{
                font-size: 6px;
            }}
            .summary-table {{
                font-size: 12px;
            }}
            .summary-table th,
            .summary-table td {{
                padding: 6px;
            }}
        }}
    </style>
</head>
<body>
<div class="container">
    <div class="header">
        <div class="ascii-container">
            <div class="ascii-art">{self.ascii_art}</div>
        </div>
        <div class="quote-container" id="quoteContainer">
            <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
            <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
        </div>
    </div>

    <div class="report-section">
        <h2>📊 Batch ABRicate Summary - {database.upper()} Database</h2>
        <div class="stats-grid">
            <div class="stat-card"><div class="stat-value">{len(hits)}</div><div class="stat-label">TOTAL HITS</div></div>
            <div class="stat-card"><div class="stat-value">{len(unique_genomes)}</div><div class="stat-label">GENOMES ANALYZED</div></div>
            <div class="stat-card"><div class="stat-value">{len(set(hit['gene'] for hit in hits))}</div><div class="stat-label">UNIQUE GENES</div></div>
            <div class="stat-card"><div class="stat-value">{datetime.now().strftime('%Y-%m-%d')}</div><div class="stat-label">ANALYSIS DATE</div></div>
        </div>
    </div>

    <div class="report-section">
        <h2>🔍 Genes by Genome</h2>
        <div class="table-responsive">
            <table class="summary-table">
                <thead><tr><th>Genome</th><th>Gene Count</th><th>Genes Detected</th></tr></thead>
                <tbody>
"""
        for genome in unique_genomes:
            genes = genes_per_genome.get(genome, set())
            gene_list = ", ".join(sorted(genes))
            html += f'<tr class="present"><td><strong>{genome}</strong></td><td>{len(genes)}</td><td>{gene_list}</td></tr>\n'

        html += """
                </tbody>
            </table>
        </div>
    </div>

    <div class="report-section">
        <h2>📈 Gene Frequency</h2>
        <div class="table-responsive">
            <table class="summary-table">
                <thead><tr><th>Gene</th><th>Frequency</th><th>Genomes</th></tr></thead>
                <tbody>
"""
        for gene, freq in sorted(gene_freq.items(), key=lambda x: x[1], reverse=True):
            # Find which genomes carry this gene
            genomes_with_gene = sorted(set(hit['genome'] for hit in hits if hit['gene'] == gene))
            genome_list = ", ".join(genomes_with_gene)
            html += f'<tr><td><strong>{gene}</strong></td><td>{freq}</td><td>{genome_list}</td></tr>\n'

        html += f"""
                </tbody>
            </table>
        </div>
    </div>

    <div class="footer">
        <p><strong>EnteroScope</strong> - Batch ABRicate Analysis Summary</p>
        <p class="timestamp">Generated: {current_time}</p>
        <div class="authorship">
            <p><strong>Author:</strong> Brown Beckley | <strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub/enteroscope" target="_blank">https://github.com/bbeckley-hub/enteroscope</a></p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
        </div>
    </div>
</div>

<script>
    const quotes = {json.dumps(self.science_quotes)};
    const quoteContainer = document.getElementById('quoteContainer');
    const quoteText = document.getElementById('quoteText');
    const quoteAuthor = document.getElementById('quoteAuthor');
    let idx = 0;
    function rotate() {{
        quoteContainer.style.opacity = '0';
        setTimeout(() => {{
            const q = quotes[idx];
            quoteText.innerHTML = '"' + q.text + '"';
            quoteAuthor.innerHTML = '— ' + q.author;
            quoteContainer.style.opacity = '1';
            idx = (idx + 1) % quotes.length;
        }}, 500);
    }}
    setInterval(rotate, 10000);
    document.addEventListener('DOMContentLoaded', rotate);
</script>
</body>
</html>"""
        out_file = os.path.join(out_base, f"enteroscope_{database}_summary_report.html")
        with open(out_file, 'w') as f:
            f.write(html)
        self.logger.info(f"Batch HTML summary: {out_file}")

    # ---------- Master JSON summary ----------
    def create_master_json_summary(self, all_results: Dict, out_base: str):
        self.logger.info("Creating master JSON summary...")
        master = {
            "metadata": self.metadata,
            "total_genomes": len(all_results),
            "databases": self.required_databases,
            "genome_summary": {},
            "overall": {"total_hits": 0, "carbapenemase_positive": 0, "esbl_positive": 0, "colistin_positive": 0}
        }
        all_hits_master = []
        for gname, gres in all_results.items():
            hits = [h for db_res in gres['results'].values() for h in db_res['hits']]
            all_hits_master.extend(hits)
            analysis = self.analyze_enterobacter_resistance(hits)
            master["genome_summary"][gname] = {
                "total_hits": len(hits),
                "carbapenemase_positive": analysis['carbapenemase_status'] == 'positive',
                "esbl_positive": analysis['esbl_status'] == 'positive',
                "colistin_positive": analysis['colistin_resistance'] == 'positive',
                "critical_carbapenemase_genes": [g['gene'] for g in analysis['critical_carbapenemase_genes']],
                "critical_esbl_genes": [g['gene'] for g in analysis['critical_esbl_genes']],
                "critical_colistin_genes": [g['gene'] for g in analysis['critical_colistin_genes']],
                "high_risk_resistance_genes": [g['gene'] for g in analysis['high_risk_resistance_genes']]
            }
            if analysis['carbapenemase_status'] == 'positive': master["overall"]["carbapenemase_positive"] += 1
            if analysis['esbl_status'] == 'positive': master["overall"]["esbl_positive"] += 1
            if analysis['colistin_resistance'] == 'positive': master["overall"]["colistin_positive"] += 1
        master["overall"]["total_hits"] = len(all_hits_master)
        out_file = os.path.join(out_base, "enteroscope_abricate_master_summary.json")
        with open(out_file, 'w') as f:
            json.dump(master, f, indent=2)
        self.logger.info(f"Master JSON: {out_file}")

    # ---------- Single and multiple genome processing ----------
    def process_single_genome(self, genome_file: str, out_base: str = "enteroscope_abricate_results") -> Dict:
        gname = Path(genome_file).stem
        out_dir = os.path.join(out_base, gname)
        os.makedirs(out_dir, exist_ok=True)
        self.logger.info(f"Processing {gname}")
        results = {}
        for db in self.required_databases:
            results[db] = self.run_abricate_single_db(genome_file, db, out_dir)
        #self.create_comprehensive_html_report(gname, results, out_dir)
        return {'genome': gname, 'results': results, 'total_hits': sum(r['hit_count'] for r in results.values())}

    def process_multiple_genomes(self, pattern: str, out_base: str = "enteroscope_abricate_results") -> Dict:
        if not self.check_abricate_installed():
            raise RuntimeError("ABRicate not installed")
        self.setup_abricate_databases()
        fasta_files = []
        for ext in ['*.fna', '*.fasta', '*.fa', '*.faa', '*.fn']:
            fasta_files.extend(glob.glob(pattern) if '*' in pattern else glob.glob(pattern + ext))
        fasta_files = list(set(fasta_files))
        if not fasta_files:
            raise FileNotFoundError(f"No FASTA files matching '{pattern}'")
        self.logger.info(f"Found {len(fasta_files)} genomes: {[Path(f).name for f in fasta_files]}")
        os.makedirs(out_base, exist_ok=True)
        all_results = {}
        if len(fasta_files) > 1 and self.cpus > 1:
            with ThreadPoolExecutor(max_workers=self.cpus) as ex:
                futures = {ex.submit(self.process_single_genome, f, out_base): f for f in fasta_files}
                for fut in as_completed(futures):
                    try:
                        res = fut.result()
                        all_results[res['genome']] = res
                    except Exception as e:
                        self.logger.error(f"Failed {futures[fut]}: {e}")
        else:
            for f in fasta_files:
                try:
                    res = self.process_single_genome(f, out_base)
                    all_results[res['genome']] = res
                except Exception as e:
                    self.logger.error(f"Failed {f}: {e}")
        self.create_database_summaries(all_results, out_base)
        self.create_master_json_summary(all_results, out_base)
        return all_results


def main():
    parser = argparse.ArgumentParser(
        description='EnteroScope ABRicate - Maximum Speed (capped at 64 cores)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python enteroscope_abricate.py "*.fna"
  python enteroscope_abricate.py "*.fasta" --cpus 64 --output my_results
  python enteroscope_abricate.py "GCF_*.fna" --cpus 32
        """
    )
    parser.add_argument('pattern', help='Glob pattern for FASTA files (e.g., "*.fna")')
    parser.add_argument('--cpus', '-c', type=int, default=None, help='Number of CPU cores (max 64)')
    parser.add_argument('--output', '-o', default='enteroscope_abricate_results', help='Output directory')
    args = parser.parse_args()

    print("\n" + "="*80)
    print("""
███████╗███╗   ██╗████████╗███████╗██████╗  ██████╗ ███████╗ ██████╗  ██████╗ ██████╗ ███████╗
██╔════╝████╗  ██║╚══██╔══╝██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔════╝ ██╔═══██╗██╔══██╗██╔════╝
█████╗  ██╔██╗ ██║   ██║   █████╗  ██████╔╝██║   ██║███████╗██║      ██║   ██║██████╔╝█████╗  
██╔══╝  ██║╚██╗██║   ██║   ██╔══╝  ██╔══██╗██║   ██║╚════██║██║      ██║   ██║██╔═══╝ ██╔══╝  
███████╗██║ ╚████║   ██║   ███████╗██║  ██║╚██████╔╝███████║╚██████╗ ╚██████╔╝██║     ███████╗
╚══════╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝ ╚═════╝  ╚═════╝ ╚═╝     ╚══════╝
""")
    print("EnteroScope ABRicate - Maximum Speed Mode (capped at 64 cores)")
    print(f"Author: Brown Beckley | GitHub: https://github.com/bbeckley-hub/enteroscope")
    print("="*80 + "\n")

    executor = EnteroAbricateExecutor(cpus=args.cpus)
    try:
        results = executor.process_multiple_genomes(args.pattern, args.output)
        executor.logger.info(f"\n✅ Analysis complete. Results in: {args.output}")
        executor.logger.info(f"   Processed {len(results)} genomes.")
    except Exception as e:
        executor.logger.error(f"Fatal error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()