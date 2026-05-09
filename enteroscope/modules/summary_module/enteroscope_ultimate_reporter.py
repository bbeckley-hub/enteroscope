#!/usr/bin/env python3
"""
EnteroScope ULTIMATE REPORTER - Enterobacter cloacae Complex
Advanced HTML Parser with Gene-Centric Cross-Genome Analysis
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
GitHub: https://github.com/bbeckley-hub/enteroscope
Version: 1.0.0
Date: 2026
"""

import os
import sys
import json
import re
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any, Optional
from datetime import datetime
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

from bs4 import BeautifulSoup

# =============================================================================
# HTML PARSER CLASS 
# =============================================================================
class EnteroUltimateHTMLParser:
    """Ultimate HTML parser for EnteroScope reports (Enterobacter cloacae complex)"""
    
    def __init__(self):
        self.all_databases = [
            'amrfinder', 'card', 'resfinder', 'argannot', 'megares', 'bacmet2',
            'vfdb', 'ecoli_vf', 'plasmidfinder', 'ecoh', 'ncbi'
        ]
        self.db_name_mapping = {
            'enteroscope_amrfinder': 'amrfinder',
            'enteroscope_card': 'card',
            'enteroscope_resfinder': 'resfinder',
            'enteroscope_argannot': 'argannot',
            'enteroscope_megares': 'megares',
            'enteroscope_bacmet2': 'bacmet2',
            'enteroscope_vfdb': 'vfdb',
            'enteroscope_ecoli_vf': 'ecoli_vf',
            'enteroscope_plasmidfinder': 'plasmidfinder',
            'enteroscope_ecoh': 'ecoh',
            'enteroscope_ncbi': 'ncbi'
        }
    
    def normalize_sample_id(self, sample_id: str) -> str:
        sample = str(sample_id)
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        extensions = ['.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', '.txt', '.tsv', '.csv']
        for ext in extensions:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
                break
        # Remove trailing .fna again (some have double)
        for ext in extensions:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
                break
        return sample.strip()
    
    def parse_html_table(self, html_content: str, table_index: int = 0) -> pd.DataFrame:
        try:
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if not tables or table_index >= len(tables):
                return pd.DataFrame()
            table = tables[table_index]
            rows = table.find_all('tr')
            if not rows:
                return pd.DataFrame()
            headers = [th.get_text().strip() for th in rows[0].find_all(['th', 'td'])]
            data = []
            for row in rows[1:]:
                cols = row.find_all(['td', 'th'])
                if cols:
                    row_data = [col.get_text().strip() for col in cols]
                    while len(row_data) < len(headers):
                        row_data.append('')
                    if len(row_data) > len(headers):
                        row_data = row_data[:len(headers)]
                    data.append(row_data)
            if not data:
                return pd.DataFrame()
            df = pd.DataFrame(data, columns=headers)
            df.columns = [col.strip().replace('\n', ' ') for col in df.columns]
            return df
        except Exception as e:
            print(f"  ⚠️ Table parsing error: {e}")
            return pd.DataFrame()
    
    def parse_mlst_report(self, file_path: Path) -> Dict[str, Dict]:
        print(f"  🧬 Parsing MLST: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            df = self.parse_html_table(html_content, 0)
            if df.empty:
                return {}
            sample_col = None
            for col in df.columns:
                if 'sample' in col.lower() or 'genome' in col.lower() or 'filename' in col.lower():
                    sample_col = col
                    break
            if not sample_col:
                sample_col = df.columns[0]
            st_col = None
            for col in df.columns:
                if col.upper() == 'ST' or col.lower() == 'st':
                    st_col = col
                    break
            if st_col is None:
                for col in df.columns:
                    if 'st' in col.lower() and 'status' not in col.lower():
                        st_col = col
                        break
            df['normalized_sample'] = df[sample_col].apply(self.normalize_sample_id)
            results = {}
            for _, row in df.iterrows():
                sample = row['normalized_sample']
                if pd.isna(sample) or not sample:
                    continue
                st = 'ND'
                if st_col is not None and pd.notna(row.get(st_col)):
                    st_val = str(row[st_col]).strip()
                    if st_val and st_val.upper() not in ['', 'NAN', 'NONE', 'ND', 'UNKNOWN']:
                        match = re.search(r'ST(\d+)', st_val, re.IGNORECASE)
                        if match:
                            st = match.group(1)
                        else:
                            num_match = re.search(r'\d+', st_val)
                            if num_match:
                                st = num_match.group()
                allele_profile = 'ND'
                for col in df.columns:
                    if 'allele' in col.lower() and 'profile' in col.lower():
                        if pd.notna(row.get(col)):
                            allele_val = str(row[col]).strip()
                            if allele_val and allele_val.lower() not in ['', 'nan', 'none', 'nd']:
                                allele_profile = allele_val
                        break
                results[sample] = {'ST': st, 'Allele_Profile': allele_profile}
            print(f"    ✓ Found {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing MLST: {e}")
            return {}
    
    def parse_qc_report(self, file_path: Path) -> Dict[str, Dict]:
        print(f"  🧬 Parsing FASTA QC: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            df = self.parse_html_table(html, 0)
            if df.empty:
                return {}
            sample_col = None
            for col in df.columns:
                if 'filename' in col.lower() or 'sample' in col.lower() or col == df.columns[0]:
                    sample_col = col
                    break
            if not sample_col:
                return {}
            results = {}
            for _, row in df.iterrows():
                sample_raw = row[sample_col]
                if not sample_raw:
                    continue
                sample = self.normalize_sample_id(sample_raw)
                qc_data = {}
                for col in df.columns:
                    if col == sample_col:
                        continue
                    val = row[col]
                    if pd.isna(val) or val == '' or val == 'ND':
                        qc_data[col] = 'ND'
                    else:
                        cleaned = str(val).replace('%', '').replace(',', '').strip()
                        try:
                            qc_data[col] = float(cleaned)
                        except:
                            qc_data[col] = str(val)
                results[sample] = qc_data
            print(f"    ✓ Parsed {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing QC: {e}")
            return {}
    
    def parse_amrfinder_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        print(f"  🧬 Parsing AMRfinder: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return {}, {}
            gene_frequencies = {}
            df_freq = self.parse_html_table(str(tables[1]), 0)
            if not df_freq.empty and 'Gene' in df_freq.columns:
                for _, row in df_freq.iterrows():
                    gene = str(row['Gene']).strip()
                    if not gene:
                        continue
                    count = 0
                    freq_str = str(row.get('Frequency', '0')).strip()
                    match = re.search(r'(\d+)', freq_str)
                    if match:
                        count = int(match.group(1))
                    percentage = (count / total_samples) * 100 if total_samples > 0 else 0
                    genomes = []
                    if 'Genomes' in df_freq.columns and pd.notna(row.get('Genomes')):
                        genomes_str = str(row['Genomes'])
                        genomes = [self.normalize_sample_id(g.strip()) for g in genomes_str.split(',') if g.strip()]
                    gene_frequencies[gene] = {
                        'count': count,
                        'percentage': round(percentage, 2),
                        'frequency_display': f"{count} ({percentage:.1f}%)",
                        'genomes': genomes,
                        'database': 'amrfinder'
                    }
            genes_by_genome = {}
            df_genomes = self.parse_html_table(str(tables[0]), 0)
            if not df_genomes.empty and 'Genome' in df_genomes.columns:
                for _, row in df_genomes.iterrows():
                    sample = self.normalize_sample_id(row['Genome'])
                    if not sample:
                        continue
                    genes = []
                    if pd.notna(row.get('Genes Detected')):
                        gene_str = str(row['Genes Detected'])
                        genes = [g.strip() for g in gene_str.split() if g.strip() and not g.startswith('Showing') and not '(' in g and not ')' in g]
                    genes_by_genome[sample] = genes
            print(f"    ✓ Found {len(genes_by_genome)} samples, {len(gene_frequencies)} genes")
            return genes_by_genome, gene_frequencies
        except Exception as e:
            print(f"    ❌ Error parsing AMRfinder: {e}")
            return {}, {}
    
    def parse_abricate_database_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        print(f"  🧬 Parsing ABRicate: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return {}, {}
            db_name = 'unknown'
            filename = str(file_path.name).lower()
            for key, value in self.db_name_mapping.items():
                if key in filename:
                    db_name = value
                    break
            gene_frequencies = {}
            df_freq = self.parse_html_table(str(tables[1]), 0)
            if not df_freq.empty and 'Gene' in df_freq.columns:
                for _, row in df_freq.iterrows():
                    gene_full = str(row['Gene']).strip()
                    if not gene_full:
                        continue
                    gene = re.sub(r'^\([^)]+\)', '', gene_full).strip()
                    if not gene:
                        gene = gene_full
                    count = 0
                    freq_str = str(row.get('Frequency', '0')).strip()
                    match = re.search(r'(\d+)', freq_str)
                    if match:
                        count = int(match.group(1))
                    percentage = (count / total_samples) * 100 if total_samples > 0 else 0
                    genomes = []
                    if 'Genomes' in df_freq.columns and pd.notna(row.get('Genomes')):
                        genomes_str = str(row['Genomes'])
                        genomes = [self.normalize_sample_id(g.strip()) for g in genomes_str.split(',') if g.strip()]
                    gene_frequencies[gene] = {
                        'count': count,
                        'percentage': round(percentage, 2),
                        'frequency_display': f"{count} ({percentage:.1f}%)",
                        'genomes': genomes,
                        'database': db_name,
                        'full_name': gene_full
                    }
            genes_by_genome = {}
            df_genomes = self.parse_html_table(str(tables[0]), 0)
            if not df_genomes.empty:
                sample_col = None
                for col in df_genomes.columns:
                    if any(keyword in col.lower() for keyword in ['genome', 'sample', 'id']):
                        sample_col = col
                        break
                if sample_col:
                    for _, row in df_genomes.iterrows():
                        sample = self.normalize_sample_id(row[sample_col])
                        if not sample:
                            continue
                        genes = []
                        genes_col = None
                        for col in df_genomes.columns:
                            if any(keyword in col.lower() for keyword in ['genes', 'detected']):
                                genes_col = col
                                break
                        if genes_col and pd.notna(row.get(genes_col)):
                            gene_str = str(row[genes_col])
                            genes = [re.sub(r'^\([^)]+\)', '', g).strip() for g in gene_str.split(',') if g.strip()]
                        genes_by_genome[sample] = genes
            print(f"    ✓ {db_name.upper()}: {len(genes_by_genome)} samples, {len(gene_frequencies)} genes")
            return genes_by_genome, gene_frequencies
        except Exception as e:
            print(f"    ❌ Error parsing ABRicate report: {e}")
            return {}, {}
    
    def parse_plasmidfinder_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        print(f"  🧬 Parsing PlasmidFinder: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return {}, {}
            plasmid_frequencies = {}
            df_freq = self.parse_html_table(str(tables[1]), 0)
            if not df_freq.empty and 'Gene' in df_freq.columns:
                for _, row in df_freq.iterrows():
                    gene_full = str(row['Gene']).strip()
                    if not gene_full:
                        continue
                    gene = self._clean_plasmid_gene_name(gene_full)
                    count = 0
                    freq_str = str(row.get('Frequency', '0')).strip()
                    match = re.search(r'(\d+)', freq_str)
                    if match:
                        count = int(match.group(1))
                    percentage = (count / total_samples) * 100 if total_samples > 0 else 0
                    genomes = []
                    if 'Genomes' in df_freq.columns and pd.notna(row.get('Genomes')):
                        genomes_str = str(row['Genomes'])
                        genomes = [self.normalize_sample_id(g.strip()) for g in genomes_str.split(',') if g.strip()]
                    plasmid_type = self._categorize_plasmid(gene)
                    plasmid_frequencies[gene] = {
                        'count': count,
                        'percentage': round(percentage, 2),
                        'frequency_display': f"{count} ({percentage:.1f}%)",
                        'genomes': genomes,
                        'full_name': gene_full,
                        'plasmid_type': plasmid_type,
                        'database': 'plasmidfinder'
                    }
            plasmids_by_genome = {}
            df_genomes = self.parse_html_table(str(tables[0]), 0)
            if not df_genomes.empty:
                sample_col = None
                for col in df_genomes.columns:
                    if any(keyword in col.lower() for keyword in ['genome', 'sample', 'id']):
                        sample_col = col
                        break
                if sample_col:
                    for _, row in df_genomes.iterrows():
                        sample = self.normalize_sample_id(row[sample_col])
                        if not sample:
                            continue
                        plasmids = []
                        genes_col = None
                        for col in df_genomes.columns:
                            if any(keyword in col.lower() for keyword in ['genes', 'detected']):
                                genes_col = col
                                break
                        if genes_col and pd.notna(row.get(genes_col)):
                            gene_str = str(row[genes_col])
                            plasmids = [self._clean_plasmid_gene_name(g.strip()) for g in gene_str.split(',') if g.strip()]
                        plasmids_by_genome[sample] = plasmids
            print(f"    ✓ PlasmidFinder: {len(plasmids_by_genome)} samples, {len(plasmid_frequencies)} plasmid markers")
            return plasmids_by_genome, plasmid_frequencies
        except Exception as e:
            print(f"    ❌ Error parsing PlasmidFinder: {e}")
            return {}, {}
    
    def _clean_plasmid_gene_name(self, gene_name: str) -> str:
        gene = re.sub(r'_(\d+)$', '', gene_name)
        plasmid_match = re.search(r'\((.*?)\)', gene)
        if plasmid_match:
            plasmid_name = plasmid_match.group(1)
            base_gene = re.sub(r'\(.*?\)', '', gene).strip()
            if base_gene:
                gene = f"{base_gene}({plasmid_name})"
            else:
                gene = plasmid_name
        return gene.strip()
    
    def _categorize_plasmid(self, gene_name: str) -> str:
        gene_lower = gene_name.lower()
        if any(rep in gene_lower for rep in ['rep', 'inc', 'rep_']):
            if 'coli' in gene_lower or 'col' in gene_lower:
                return 'Colicin plasmid'
            elif 'broad' in gene_lower or 'multihost' in gene_lower:
                return 'Broad-host-range plasmid'
            else:
                return 'Replication protein'
        elif any(mob in gene_lower for mob in ['mob', 'tra', 'conj']):
            return 'Mobility/conjugation'
        elif 'col' in gene_lower and 'coli' not in gene_lower:
            return 'Colicin plasmid'
        elif 'inc' in gene_lower and gene_lower.startswith('inc'):
            inc_match = re.search(r'inc([a-z0-9]+)', gene_lower)
            if inc_match:
                inc_type = inc_match.group(1).upper()
                return f'Incompatibility group {inc_type}'
            return 'Incompatibility group'
        else:
            return 'Other plasmid'


# =============================================================================
# DATA ANALYZER CLASS (Enterobacter-specific gene sets)
# =============================================================================
class EnteroDataAnalyzer:
    def __init__(self):
        # Critical carbapenemases
        self.critical_carbapenemases = {
            'blaKPC', 'blaNDM', 'blaVIM', 'blaIMP', 'blaOXA-48', 'blaOXA-181', 'blaOXA-232',
            'blaGES', 'blaIMI', 'blaSME', 'blaNMC', 'blaCcrA', 'blaCphA',
            'KPC', 'NDM', 'VIM', 'IMP', 'OXA-48', 'OXA-181', 'OXA-232',
            'GES', 'IMI', 'SME', 'NMC', 'CcrA', 'CphA'
        }
        # ESBLs / AmpC
        self.critical_esbls = {
            'blaCTX-M', 'blaSHV', 'blaTEM', 'blaPER', 'blaVEB', 'blaGES', 'blaBEL',
            'CTX-M', 'SHV', 'TEM', 'PER', 'VEB', 'GES', 'BEL', 'BES',
            'blaCMY', 'blaDHA', 'blaFOX', 'blaMOX', 'blaACC', 'blaACT', 'blaMIR',
            'CMY', 'DHA', 'FOX', 'MOX', 'ACC', 'ACT', 'MIR'
        }
        # Colistin
        self.critical_colistin = {
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10',
            'mcr1', 'mcr2', 'mcr3', 'mcr4', 'mcr5', 'mcr6', 'mcr7', 'mcr8', 'mcr9', 'mcr10',
            'pmrA', 'pmrB', 'pmrC', 'lpxA', 'lpxC', 'lpxD', 'arnA', 'arnB', 'arnC', 'arnD', 'eptA'
        }
        # Tigecycline
        self.critical_tigecycline = {
            'tet(X)', 'tet(X1)', 'tet(X2)', 'tet(X3)', 'tet(X4)', 'tet(X5)', 'tet(X6)',
            'tet(39)', 'tet(A)', 'tet(B)', 'tet(C)', 'tet(D)', 'tet(E)', 'tet(G)', 'tet(H)',
            'adeS', 'adeR', 'adeA', 'adeB', 'adeC', 'adeJ', 'adeK', 'adeN', 'adeT'
        }
        # Biofilm
        self.critical_biofilm = {
            'ompA', 'csuA', 'csuB', 'csuC', 'csuD', 'csuE', 'csuA/B', 'csuC/D/E',
            'bfmR', 'bfmS', 'abaI', 'abaR', 'pilA', 'pilB', 'pilC', 'pilD', 'pilE', 'pilF',
            'ptk', 'epsA', 'pgaA', 'pgaB', 'pgaC', 'pgaD', 'bap', 'csgA', 'csgB', 'csgD'
        }
        # Efflux pumps
        self.critical_efflux = {
            'adeA', 'adeB', 'adeC', 'adeF', 'adeG', 'adeH', 'adeI', 'adeJ', 'adeK', 'adeL', 'adeM', 'adeN',
            'abeM', 'abeS', 'acrA', 'acrB', 'acrD', 'acrF', 'tolC', 'mexA', 'mexB', 'mexC', 'mexD', 'mexE', 'mexF',
            'oprM', 'oprN', 'oprJ', 'mdtA', 'mdtB', 'mdtC', 'emrA', 'emrB', 'emrD', 'emrE', 'emrK', 'emrY',
            'mdx', 'acrAB', 'acrEF', 'mdfA', 'cmlA', 'oqxA', 'oqxB', 'oqxR'
        }
        # Bacmet2 markers (biocides, heavy metals)
        self.bacmet2_markers = {
            'qacA', 'qacB', 'qacC', 'qacD', 'qacE', 'qacF', 'qacG', 'qacH', 'qacI', 'qacJ',
            'qacEA1', 'qacG2', 'qacH2', 'cepA', 'formA', 'formB', 'formC', 'oqxA', 'oqxB',
            'czcA', 'czcB', 'czcC', 'czcD', 'czcR', 'czcS',
            'merA', 'merB', 'merC', 'merD', 'merE', 'merF', 'merG', 'merH', 'merI', 'merJ', 'merP', 'merT',
            'arsA', 'arsB', 'arsC', 'arsD', 'arsE', 'arsF', 'arsG', 'arsH', 'arsI', 'arsJ', 'arsT',
            'copA', 'copB', 'copC', 'copD', 'copE', 'copF', 'copG', 'copH', 'copI', 'copJ',
            'zntA', 'zntB', 'zntC', 'zntD', 'zntE', 'zntF', 'zntG', 'zntH', 'zntI', 'zntJ',
            'chrA', 'chrB', 'chrC', 'chrD', 'chrE', 'chrF',
            'nikA', 'nikB', 'nikC', 'nikD', 'nikE', 'nikR',
            'cadA', 'cadB', 'cadC', 'cadD',
            'silA', 'silB', 'silC', 'silD', 'silE',
            'pbrA', 'pbrB', 'pbrC', 'pbrD', 'pbrR',
            'corA', 'corC', 'corR', 'zraR', 'pitA', 'nccN', 'nreB', 'fptA', 'fecE', 'fpvA', 'znuB', 'znuC', 'frnE'
        }
    
    def categorize_gene(self, gene: str) -> str:
        import re
        gl = gene.lower()
        
        # Helper: check if any pattern exists as a whole word (case‑insensitive)
        def has_word(patterns):
            for pat in patterns:
                if re.search(r'\b' + re.escape(pat.lower()) + r'\b', gl):
                    return True
            return False
        
        # Carbapenemases (exact word match to avoid "fimI" matching "IMI")
        if has_word(['kpc', 'ndm', 'vim', 'imp', 'oxa-48', 'oxa-181', 'oxa-232',
                    'ges', 'imi', 'sme', 'nmc', 'ccra', 'cpha', 'blaoxa']):
            return 'Carbapenemases'
        
        # ESBLs / AmpC
        esbl_words = ['ctx-m', 'shv', 'tem', 'per', 'veb', 'bel', 'bes',
                    'cmy', 'dha', 'fox', 'mox', 'acc', 'act', 'mir']
        if has_word(esbl_words):
            return 'ESBLs'
        
        # Colistin resistance
        col_words = ['mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6',
                    'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10', 'pmra', 'pmrb',
                    'pmrc', 'lpxa', 'lpxc', 'lpxd', 'arna', 'arnb', 'arnc',
                    'arnd', 'epta']
        if has_word(col_words):
            return 'Colistin Resistance'
        
        # Tigecycline resistance
        tig_words = ['tet(x', 'tet(39', 'tet(a', 'tet(b', 'tet(c', 'tet(d',
                    'tet(e', 'tet(g', 'tet(h', 'ades', 'ader', 'adea',
                    'adeb', 'adec', 'adej', 'adek', 'aden', 'adet']
        if has_word(tig_words):
            return 'Tigecycline Resistance'
        
        # Biofilm
        biofilm_words = ['ompa', 'csu', 'bfmr', 'bfms', 'abai', 'abar', 'pila',
                        'pilb', 'pilc', 'pild', 'pile', 'pilf', 'ptk', 'epsa',
                        'pgaa', 'pgab', 'pgac', 'pgad', 'bap', 'csga', 'csgb',
                        'csgd']
        if has_word(biofilm_words):
            return 'Biofilm Formation'
        
        # Efflux pumps
        efflux_words = ['ade', 'abem', 'abes', 'acra', 'acrb', 'acrd', 'acrf',
                        'tolc', 'mex', 'oprm', 'oprn', 'oprj', 'mdta', 'mdtb',
                        'mdtc', 'emra', 'emrb', 'emrd', 'emre', 'emrk', 'emry',
                        'mdx', 'mdfa', 'cmla', 'oqx']
        if has_word(efflux_words):
            return 'Efflux Pumps'
        
        # Bacmet (biocide/metal)
        bacmet_words = ['qac', 'cepa', 'forma', 'czc', 'mer', 'ars', 'cop',
                        'znt', 'chr', 'nik', 'cad', 'sil', 'pbr', 'cor', 'zrar',
                        'pita', 'nccn', 'nreb', 'fpta', 'fece', 'fpva', 'znub',
                        'znuc', 'frne']
        if has_word(bacmet_words):
            return 'Bacmet (Biocide/Metal)'
        
        return 'Other'
    
    def create_gene_centric_tables(self, integrated_data: Dict[str, Any], total_samples: int) -> Dict[str, Any]:
        gene_centric = {
            'amr_databases': {},
            'virulence_databases': {},
            'bacmet_databases': {},
            'plasmid_databases': {},
            'combined_gene_frequencies': [],
            'gene_categories': defaultdict(list),
            'database_stats': {}
        }
        if 'amrfinder' in integrated_data.get('gene_frequencies', {}):
            amr_data = integrated_data['gene_frequencies']['amrfinder']
            gene_list = []
            for gene, data in amr_data.items():
                category = self.categorize_gene(gene)
                gene_list.append({
                    'gene': gene, 'category': category, 'database': 'AMRfinder',
                    'count': data['count'], 'percentage': data['percentage'],
                    'frequency_display': data['frequency_display'],
                    'genomes': data['genomes']
                })
            if gene_list:
                gene_list.sort(key=lambda x: x['count'], reverse=True)
                gene_centric['amr_databases']['amrfinder'] = gene_list
                for gd in gene_list:
                    gene_centric['gene_categories'][gd['category']].append(gd)
        if 'abricate' in integrated_data.get('gene_frequencies', {}):
            abricate_data = integrated_data['gene_frequencies']['abricate']
            for db_name, db_genes in abricate_data.items():
                gene_list = []
                for gene, data in db_genes.items():
                    category = self.categorize_gene(gene)
                    gene_list.append({
                        'gene': gene, 'category': category, 'database': db_name.upper(),
                        'count': data['count'], 'percentage': data['percentage'],
                        'frequency_display': data['frequency_display'],
                        'genomes': data['genomes'], 'full_name': data.get('full_name', gene)
                    })
                if gene_list:
                    gene_list.sort(key=lambda x: x['count'], reverse=True)
                    if db_name in ['vfdb', 'ecoli_vf']:
                        gene_centric['virulence_databases'][db_name] = gene_list
                    elif db_name == 'bacmet2':
                        gene_centric['bacmet_databases'][db_name] = gene_list
                    elif db_name in ['plasmidfinder', 'ecoh']:
                        gene_centric['plasmid_databases'][db_name] = gene_list
                    else:
                        gene_centric['amr_databases'][db_name] = gene_list
                    for gd in gene_list:
                        gene_centric['gene_categories'][gd['category']].append(gd)
        # Combined
        all_genes = []
        for db_type in ['amr_databases', 'virulence_databases', 'bacmet_databases', 'plasmid_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                all_genes.extend(genes)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        gene_centric['combined_gene_frequencies'] = all_genes
        for cat in gene_centric['gene_categories']:
            gene_centric['gene_categories'][cat].sort(key=lambda x: x['count'], reverse=True)
        # Stats
        for db_type in ['amr_databases', 'virulence_databases', 'bacmet_databases', 'plasmid_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                gene_centric['database_stats'][db_name] = {
                    'total_genes': len(genes),
                    'total_occurrences': sum(g['count'] for g in genes),
                    'critical_genes': sum(1 for g in genes if g['category'] in ['Carbapenemases', 'ESBLs', 'Colistin Resistance', 'Tigecycline Resistance']),
                }
        return gene_centric
    
    def create_cross_genome_patterns(self, integrated_data: Dict[str, Any], total_samples: int) -> Dict[str, Any]:
        patterns = {
            'st_distribution': Counter(),
            'high_risk_isolates': []
        }
        samples_data = integrated_data.get('samples', {})
        gene_centric = integrated_data.get('gene_centric', {})
        sample_genes = defaultdict(set)
        for db_type in ['amr_databases', 'virulence_databases', 'bacmet_databases', 'plasmid_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for gene_data in genes:
                    for genome in gene_data['genomes']:
                        sample_genes[genome].add(gene_data['gene'])
        for sample, data in samples_data.items():
            st = data.get('mlst', {}).get('ST', 'ND')
            if st != 'ND':
                patterns['st_distribution'][st] += 1
            genes = list(sample_genes.get(sample, set()))
            carb = [g for g in genes if 'Carbapenemases' in self.categorize_gene(g)]
            col = [g for g in genes if 'Colistin Resistance' in self.categorize_gene(g)]
            tig = [g for g in genes if 'Tigecycline Resistance' in self.categorize_gene(g)]
            if carb and (col or tig):
                patterns['high_risk_isolates'].append({
                    'sample': sample, 'st': st,
                    'carbapenemases': carb, 'colistin_resistance': col,
                    'tigecycline_resistance': tig
                })
        return patterns
    
    def create_plasmid_analysis(self, integrated_data: Dict[str, Any], total_samples: int) -> Dict[str, Any]:
        plasmid_analysis = {'plasmid_databases': {}, 'plasmid_frequencies': [], 'sample_plasmid_profiles': defaultdict(list)}
        if 'plasmidfinder' not in integrated_data.get('gene_frequencies', {}):
            return plasmid_analysis
        plasmid_genes = integrated_data['gene_frequencies']['plasmidfinder']
        gene_list = []
        for gene, data in plasmid_genes.items():
            gene_list.append({
                'plasmid_marker': gene,
                'category': data.get('plasmid_type', 'Other'),
                'database': 'PlasmidFinder',
                'count': data['count'],
                'percentage': data['percentage'],
                'frequency_display': data['frequency_display'],
                'genomes': data['genomes']
            })
        if gene_list:
            gene_list.sort(key=lambda x: x['count'], reverse=True)
            plasmid_analysis['plasmid_databases']['plasmidfinder'] = gene_list
            plasmid_analysis['plasmid_frequencies'] = gene_list
        samples_data = integrated_data.get('samples', {})
        for sample in samples_data:
            if 'plasmidfinder' in integrated_data.get('gene_frequencies', {}):
                db_genes = integrated_data['gene_frequencies']['plasmidfinder']
                sample_plasmids = []
                for gene, gdata in db_genes.items():
                    if sample in gdata.get('genomes', []):
                        sample_plasmids.append(gene)
                if sample_plasmids:
                    plasmid_analysis['sample_plasmid_profiles'][sample] = sample_plasmids
        return plasmid_analysis


# =============================================================================
# HTML GENERATOR CLASS (template style, teal/cyan theme)
# =============================================================================
class EnteroHTMLGenerator:
    def __init__(self, data_analyzer: EnteroDataAnalyzer):
        self.data_analyzer = data_analyzer
        # 20 educational biology facts about Enterobacter cloacae complex
        self.educational_facts = [
            "Enterobacter cloacae complex is a group of opportunistic pathogens commonly found in hospital environments, soil, and water.",
            "Members of the Enterobacter cloacae complex are intrinsically resistant to ampicillin and amoxicillin due to a chromosomal AmpC beta-lactamase.",
            "Carbapenem-resistant Enterobacter cloacae (CREC) is an emerging threat with mortality rates exceeding 40% in vulnerable patients.",
            "The OXA-48 carbapenemase is frequently detected in Enterobacter cloacae, often carried on epidemic IncL/M plasmids that facilitate horizontal spread.",
            "Enterobacter cloacae can acquire colistin resistance via plasmid-borne mcr genes or through chromosomal mutations in the pmrAB two-component system.",
            "Biofilm formation by Enterobacter cloacae on medical devices (catheters, ventilators) is a major cause of persistent healthcare-associated infections.",
            "Tigecycline resistance in Enterobacter cloacae is often mediated by tet(X) variants or overexpression of efflux pumps like AdeABC and OqxAB.",
            "Enterobacter cloacae complex comprises at least 12 distinct species, including E. cloacae, E. hormaechei, E. asburiae, and E. kobei, with varying clinical relevance.",
            "The CRISPR-Cas system is present in some Enterobacter cloacae strains, providing adaptive immunity against bacteriophages and mobile genetic elements.",
            "Enterobacter cloacae is a leading cause of neonatal sepsis in low- and middle-income countries, often associated with contaminated intravenous fluids.",
            "The Enterobacter cloacae complex is intrinsically resistant to cephamycins (cefoxitin, cefotetan) due to inducible AmpC production.",
            "Plasmid-mediated quinolone resistance (PMQR) genes such as qnrB and aac(6')-Ib-cr are increasingly reported in Enterobacter cloacae clinical isolates.",
            "Enterobacter cloacae produces a variety of siderophores (enterobactin, aerobactin) for iron acquisition, enhancing its virulence in iron-limited host environments.",
            "Some Enterobacter cloacae strains harbor the blaIMI gene, a carbapenemase that is chromosomally encoded but can be mobilized by insertion sequences.",
            "The Enterobacter cloacae complex is a frequent colonizer of the human gut, serving as a reservoir for resistance genes that can transfer to other Enterobacteriaceae.",
            "Hospital outbreaks of Enterobacter cloacae are often associated with contaminated handwashing sinks, endoscopes, and ultrasound gel.",
            "Efflux pumps such as AcrAB-TolC and OqxAB contribute significantly to multidrug resistance in Enterobacter cloacae, including reduced susceptibility to tigecycline and chloramphenicol.",
            "The lipopolysaccharide (LPS) of Enterobacter cloacae can be modified by the arn operon, conferring resistance to polymyxins (colistin).",
            "Enterobacter cloacae possesses a type VI secretion system (T6SS) that mediates interbacterial competition and may contribute to colonization resistance.",
            "Whole-genome sequencing and MLST have revealed that sequence type ST171, ST78, ST114, and ST418 are globally disseminated clones of carbapenem-resistant Enterobacter cloacae."
        ]
    
    def get_random_fact(self):
        import random
        return random.choice(self.educational_facts)
    
    def generate_main_report(self, integrated_data: Dict[str, Any], output_dir: Path) -> str:
        print("\n🎨 Generating ULTIMATE HTML report for Enterobacter cloacae complex...")
        samples_data = integrated_data.get('samples', {})
        patterns = integrated_data.get('patterns', {})
        gene_centric = integrated_data.get('gene_centric', {})
        metadata = integrated_data.get('metadata', {})
        plasmid_analysis = integrated_data.get('plasmid_analysis', {})
        qc_data = integrated_data.get('qc_data', {})
        html = self._create_ultimate_html(
            metadata=metadata, samples_data=samples_data, patterns=patterns,
            gene_centric=gene_centric, plasmid_analysis=plasmid_analysis,
            qc_data=qc_data, total_samples=len(samples_data)
        )
        output_file = output_dir / "enteroscope_ultimate_report.html"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"    ✅ HTML report saved: {output_file}")
        return str(output_file)
    
    def _create_ultimate_html(self, **kwargs) -> str:
        metadata = kwargs['metadata']
        samples_data = kwargs['samples_data']
        patterns = kwargs['patterns']
        gene_centric = kwargs['gene_centric']
        plasmid_analysis = kwargs.get('plasmid_analysis', {})
        qc_data = kwargs.get('qc_data', {})
        total_samples = kwargs['total_samples']
        total_amr = sum(len(db) for db in gene_centric.get('amr_databases', {}).values())
        total_vir = sum(len(db) for db in gene_centric.get('virulence_databases', {}).values())
        total_bacmet = sum(len(db) for db in gene_centric.get('bacmet_databases', {}).values())
        total_plasmid = sum(len(db) for db in plasmid_analysis.get('plasmid_databases', {}).values())
        high_risk = len(patterns.get('high_risk_isolates', []))
        carb_count = sum(1 for db in gene_centric.get('amr_databases', {}).values() for g in db if g['category'] == 'Carbapenemases')
        
        # Teal/Cyan CSS
        css = """
        <style>
        :root { --summary-color:#4CAF50; --samples-color:#2196F3; --mlst-color:#FF9800; --amr-color:#F44336; --virulence-color:#E91E63; --bacmet-color:#795548; --combinations-color:#009688; --patterns-color:#FF5722; --plasmid-color:#2196F3; --databases-color:#607D8B; --qc-color:#17a2b8; --aiguide-color:#3F51B5; --calltoaction-color:#2E7D32; --export-color:#3F51B5; }
        *{margin:0;padding:0;box-sizing:border-box}
        body{font-family:'Segoe UI',Tahoma,Geneva,Verdana,sans-serif;line-height:1.6;color:#333;background:#f5f5f5;overflow-x:auto}
        .container{width:100%;max-width:100%;margin:0 auto;padding:20px;overflow-x:hidden}
        .main-header{background:linear-gradient(135deg,#0f766e 0%,#0d9488 100%);color:white;padding:30px;border-radius:15px;margin-bottom:30px;text-align:center}
        .main-header h1{font-size:2.8em;margin-bottom:10px;color:white}
        .metadata-bar{background:rgba(255,255,255,0.1);padding:15px;border-radius:10px;margin:20px 0;display:flex;justify-content:space-around;flex-wrap:wrap;gap:15px}
        .quote-container{background:rgba(0,0,0,0.4);backdrop-filter:blur(5px);padding:20px;border-radius:12px;margin:20px 0;text-align:center;border-left:4px solid #06b6d4;transition:opacity 0.5s}
        .quote-text{font-size:18px;font-style:italic;margin-bottom:10px;color:white}
        .dashboard-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:20px;margin-bottom:30px}
        .dashboard-card{background:white;padding:25px;border-radius:12px;box-shadow:0 5px 20px rgba(0,0,0,0.1);text-align:center;cursor:pointer;border-left:5px solid;transition:all 0.3s}
        .dashboard-card:hover{transform:translateY(-10px)}
        .card-number{font-size:3em;font-weight:bold;margin:15px 0;background:linear-gradient(90deg,#0f766e,#0d9488);-webkit-background-clip:text;-webkit-text-fill-color:transparent}
        .tab-navigation{display:flex;gap:5px;margin-bottom:20px;flex-wrap:wrap;background:white;padding:15px;border-radius:12px;position:sticky;top:10px;z-index:100}
        .tab-button{padding:12px 25px;background:#f5f5f5;border:none;border-radius:8px;cursor:pointer;font-weight:600;display:flex;align-items:center;gap:8px;transition:0.3s}
        .tab-button.active{color:white}
        .tab-button.summary.active{background:var(--summary-color)}
        .tab-button.samples.active{background:var(--samples-color)}
        .tab-button.mlst.active{background:var(--mlst-color)}
        .tab-button.amr.active{background:var(--amr-color)}
        .tab-button.virulence.active{background:var(--virulence-color)}
        .tab-button.bacmet.active{background:var(--bacmet-color)}
        .tab-button.combinations.active{background:var(--combinations-color)}
        .tab-button.patterns.active{background:var(--patterns-color)}
        .tab-button.plasmid.active{background:var(--plasmid-color)}
        .tab-button.databases.active{background:var(--databases-color)}
        .tab-button.qc.active{background:var(--qc-color)}
        .tab-button.aiguide.active{background:var(--aiguide-color)}
        .tab-button.calltoaction.active{background:var(--calltoaction-color)}
        .tab-button.export.active{background:var(--export-color)}
        .tab-content{display:none;background:white;padding:30px;border-radius:15px;margin-bottom:30px;animation:fadeIn 0.5s}
        .tab-content.active{display:block}
        @keyframes fadeIn{from{opacity:0;transform:translateY(20px)}to{opacity:1;transform:translateY(0)}}
        .section-header{color:#2c3e50;margin-bottom:25px;padding-bottom:15px;border-bottom:3px solid;font-size:1.8em;display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap}
        .master-scrollable-container{width:100%;overflow-x:auto;border:1px solid #e0e0e0;border-radius:8px;margin:20px 0}
        .data-table{width:100%;border-collapse:collapse;font-size:0.95em;box-shadow:0 2px 10px rgba(0,0,0,0.1);min-width:600px}
        .data-table th{background:#0f766e;color:white;padding:15px;text-align:left;white-space:nowrap;cursor:pointer}
        .data-table th:hover{background:#0d5c56}
        .data-table td{padding:12px 15px;border-bottom:1px solid #e0e0e0;vertical-align:top;word-break:break-word}
        .data-table tr:hover{background:#f8f9fa}
        .data-table td:first-child, .data-table th:first-child { min-width: 200px; white-space: nowrap; }
        .search-box{width:100%;padding:12px;margin-bottom:15px;border:2px solid #e0e0e0;border-radius:8px;font-size:1em}
        .search-box:focus{outline:none;border-color:#0f766e}
        .action-buttons{display:flex;gap:10px;margin:20px 0;flex-wrap:wrap}
        .action-btn{padding:10px 20px;border:none;border-radius:8px;cursor:pointer;font-weight:600;display:inline-flex;align-items:center;gap:8px;transition:0.3s}
        .action-btn:hover{transform:translateY(-2px)}
        .btn-primary{background:#0f766e;color:white}
        .btn-success{background:#28a745;color:white}
        .btn-warning{background:#ffc107;color:black}
        .btn-secondary{background:#6c757d;color:white}
        .badge{display:inline-block;padding:5px 15px;border-radius:20px;font-size:0.85em;font-weight:600;margin:2px}
        .badge-critical{background:#9C27B0;color:white}
        .badge-high{background:#F44336;color:white}
        .badge-medium{background:#FF9800;color:black}
        .alert-box{padding:20px;border-radius:10px;margin:20px 0;display:flex;align-items:center;gap:20px;border-left:5px solid}
        .alert-info{background:#d1ecf1;color:#0c5460;border-left-color:#17a2b8}
        .alert-danger{background:#f8d7da;color:#721c24;border-left-color:#dc3545}
        .genome-list{display:block;max-height:200px;overflow-y:auto;padding:5px;background:#f8f9fa;border-radius:5px}
        .genome-tag{display:inline-block;background:#e0f2f1;color:#0f766e;padding:3px 10px;border-radius:12px;font-size:0.85em;margin:2px}
        .genome-tag.highlight{background-color:#ffff99 !important; color:#000 !important; border:1px solid #ffc107}
        .footer{text-align:center;padding:30px;background:linear-gradient(135deg,#2c3e50,#34495e);color:white;border-radius:15px;margin-top:40px}
        .footer a{color:#ffc107;text-decoration:none}
        .footer a:hover{text-decoration:underline}
        .info-text{background:#f8f9fa;padding:15px;border-radius:8px;margin:15px 0;border-left:4px solid #17a2b8}
        .critical-cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:20px;margin:20px 0}
        .critical-card{background:#fff;padding:20px;border-radius:10px;box-shadow:0 2px 10px rgba(0,0,0,0.1);border-left:4px solid}
        .sort-icon{margin-left:5px;font-size:0.8em;opacity:0.6}
        </style>
        """
        js = """
        <script>
        function switchTab(tabName){
            document.querySelectorAll('.tab-content').forEach(t=>t.classList.remove('active'));
            document.querySelectorAll('.tab-button').forEach(b=>b.classList.remove('active'));
            document.getElementById(tabName+'-tab').classList.add('active');
            event.currentTarget.classList.add('active');
            window.location.hash=tabName;
        }
        function searchTable(tableId, searchId){
            let filter=document.getElementById(searchId).value.toUpperCase();
            let rows=document.getElementById(tableId).getElementsByTagName('tr');
            for(let i=1;i<rows.length;i++){
                let found=false;
                for(let cell of rows[i].getElementsByTagName('td')){
                    if(cell.textContent.toUpperCase().indexOf(filter)>-1){found=true;break;}
                }
                rows[i].style.display=found?'':'none';
            }
        }
        function highlightGenome(tableId, genomeSearchId){
            let filter=document.getElementById(genomeSearchId).value.toUpperCase().trim();
            let table=document.getElementById(tableId);
            let allTags=table.querySelectorAll('.genome-tag');
            allTags.forEach(tag=>tag.classList.remove('highlight'));
            if(filter==='') return;
            allTags.forEach(tag=>{
                if(tag.textContent.toUpperCase().indexOf(filter)>-1){
                    tag.classList.add('highlight');
                }
            });
        }
        function exportTableToCSV(tableId, filename){
            let rows=document.getElementById(tableId).querySelectorAll('tr');
            let csv=[];
            for(let row of rows){
                let cols=row.querySelectorAll('td,th');
                let rowData=Array.from(cols).map(c=>'"'+c.innerText.replace(/"/g,'""')+'"');
                csv.push(rowData.join(','));
            }
            let blob=new Blob([csv.join('\\n')],{type:'text/csv'});
            let a=document.createElement('a');
            a.href=URL.createObjectURL(blob);
            a.download=filename;
            a.click();
            URL.revokeObjectURL(a.href);
        }
        function printSection(sectionId){
            let content=document.getElementById(sectionId);
            let win=window.open('','_blank');
            win.document.write('<html><head><title>Print</title><style>'+document.querySelector('style').innerHTML+'</style></head><body>'+content.innerHTML+'</body></html>');
            win.document.close();
            win.print();
        }
        function sortTable(tableId, colIndex, type){
            let table=document.getElementById(tableId);
            let tbody=table.tBodies[0];
            let rows=Array.from(tbody.rows);
            let isAscending=table.getAttribute('data-sort-dir')!=='asc';
            rows.sort((a,b)=>{
                let aVal=a.cells[colIndex].innerText.trim();
                let bVal=b.cells[colIndex].innerText.trim();
                if(type==='number'){
                    aVal=parseFloat(aVal.replace(/,/g,''))||0;
                    bVal=parseFloat(bVal.replace(/,/g,''))||0;
                    return isAscending?aVal-bVal:bVal-aVal;
                }else{
                    return isAscending?aVal.localeCompare(bVal):bVal.localeCompare(aVal);
                }
            });
            tbody.append(...rows);
            table.setAttribute('data-sort-dir',isAscending?'asc':'desc');
            let headers=table.querySelectorAll('th');
            headers.forEach((th,idx)=>{
                let icon=th.querySelector('.sort-icon');
                if(icon) icon.innerHTML='⇅';
            });
            let currentHeader=headers[colIndex];
            let icon=currentHeader.querySelector('.sort-icon');
            if(icon) icon.innerHTML=isAscending?'↑':'↓';
        }
        document.addEventListener('DOMContentLoaded',function(){
            let hash=window.location.hash.substring(1);
            if(hash){
                let btn=document.querySelector(`.tab-button.${hash}`);
                if(btn) btn.click();
            } else document.querySelector('.tab-button').click();
            document.querySelectorAll('.data-table').forEach(table=>{
                let headers=table.querySelectorAll('th');
                headers.forEach((th,idx)=>{
                    let type=th.getAttribute('data-sort')||'string';
                    th.style.cursor='pointer';
                    th.addEventListener('click',()=>sortTable(table.id,idx,type));
                    let icon=document.createElement('span');
                    icon.className='sort-icon';
                    icon.innerHTML='⇅';
                    th.appendChild(icon);
                });
            });
        });
        </script>
        """
        # Initial random fact with emoji prefix
        initial_fact = self.get_random_fact()
        
        html = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>EnteroScope Ultimate Report - Enterobacter cloacae Complex</title>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
{css}{js}</head>
<body><div class="container">
<div class="main-header"><h1><i class="fas fa-bacterium"></i> EnteroScope Ultimate Analysis Report</h1>
<p>Gene-Centric Cross-Genome Analysis for <em>Enterobacter cloacae</em> Complex</p>
<div class="metadata-bar"><div class="metadata-item"><i class="fas fa-calendar"></i> {metadata.get('analysis_date','Unknown')}</div>
<div class="metadata-item"><i class="fas fa-database"></i> Samples: {total_samples}</div>
<div class="metadata-item"><i class="fas fa-vial"></i> Pathogen: Enterobacter cloacae complex</div>
<div class="metadata-item"><i class="fas fa-university"></i> University of Ghana Medical School</div></div>
<div class="quote-container" id="quoteContainer"><div class="quote-text" id="quoteText">🧬 {initial_fact}</div></div></div>
<div class="dashboard-grid">
<div class="dashboard-card card-summary" onclick="switchTab('summary')"><div class="card-number">{total_samples}</div><div class="card-label">Total Samples</div><i class="fas fa-vial fa-2x" style="color:var(--summary-color)"></i></div>
<div class="dashboard-card card-mlst" onclick="switchTab('mlst')"><div class="card-number">{len(patterns.get('st_distribution',{}))}</div><div class="card-label">Sequence Types</div><i class="fas fa-code-branch fa-2x" style="color:var(--mlst-color)"></i></div>
<div class="dashboard-card card-amr" onclick="switchTab('amr')"><div class="card-number">{total_amr}</div><div class="card-label">AMR Genes</div><i class="fas fa-biohazard fa-2x" style="color:var(--amr-color)"></i></div>
<div class="dashboard-card card-virulence" onclick="switchTab('virulence')"><div class="card-number">{total_vir}</div><div class="card-label">Virulence Genes</div><i class="fas fa-virus fa-2x" style="color:var(--virulence-color)"></i></div>
<div class="dashboard-card card-bacmet" onclick="switchTab('bacmet')"><div class="card-number">{total_bacmet}</div><div class="card-label">Bacmet</div><i class="fas fa-flask fa-2x" style="color:var(--bacmet-color)"></i></div>
<div class="dashboard-card card-plasmid" onclick="switchTab('plasmid')"><div class="card-number">{total_plasmid}</div><div class="card-label">Plasmid Markers</div><i class="fas fa-dna fa-2x" style="color:var(--plasmid-color)"></i></div>
<div class="dashboard-card card-patterns" onclick="switchTab('patterns')"><div class="card-number">{high_risk}</div><div class="card-label">High-Risk Isolates</div><i class="fas fa-exclamation-triangle fa-2x" style="color:var(--patterns-color)"></i></div>
</div>
<div class="tab-navigation">
<button class="tab-button summary active" onclick="switchTab('summary')"><i class="fas fa-chart-pie"></i> Summary</button>
<button class="tab-button samples" onclick="switchTab('samples')"><i class="fas fa-list-alt"></i> Sample Overview</button>
<button class="tab-button qc" onclick="switchTab('qc')"><i class="fas fa-chart-line"></i> FASTA QC</button>
<button class="tab-button mlst" onclick="switchTab('mlst')"><i class="fas fa-code-branch"></i> MLST</button>
<button class="tab-button amr" onclick="switchTab('amr')"><i class="fas fa-biohazard"></i> AMR</button>
<button class="tab-button virulence" onclick="switchTab('virulence')"><i class="fas fa-virus"></i> Virulence</button>
<button class="tab-button bacmet" onclick="switchTab('bacmet')"><i class="fas fa-flask"></i> Bacmet</button>
<button class="tab-button plasmid" onclick="switchTab('plasmid')"><i class="fas fa-dna"></i> Plasmids</button>
<button class="tab-button patterns" onclick="switchTab('patterns')"><i class="fas fa-project-diagram"></i> Patterns</button>
<button class="tab-button databases" onclick="switchTab('databases')"><i class="fas fa-database"></i> Database Metrics</button>
<button class="tab-button aiguide" onclick="switchTab('aiguide')"><i class="fas fa-robot"></i> AI Guide</button>
<button class="tab-button calltoaction" onclick="switchTab('calltoaction')"><i class="fas fa-globe"></i> Call to Action</button>
<button class="tab-button export" onclick="switchTab('export')"><i class="fas fa-download"></i> Export</button>
</div>
<div id="summary-tab" class="tab-content active">{self._summary_section(kwargs, carb_count, total_bacmet)}</div>
<div id="samples-tab" class="tab-content">{self._samples_section(kwargs, qc_data)}</div>
<div id="qc-tab" class="tab-content">{self._qc_section(kwargs)}</div>
<div id="mlst-tab" class="tab-content">{self._mlst_section(kwargs)}</div>
<div id="amr-tab" class="tab-content">{self._amr_section(kwargs)}</div>
<div id="virulence-tab" class="tab-content">{self._virulence_section(kwargs)}</div>
<div id="bacmet-tab" class="tab-content">{self._bacmet_section(kwargs)}</div>
<div id="plasmid-tab" class="tab-content">{self._plasmid_section(kwargs)}</div>
<div id="patterns-tab" class="tab-content">{self._patterns_section(kwargs)}</div>
<div id="databases-tab" class="tab-content">{self._databases_section(kwargs)}</div>
<div id="aiguide-tab" class="tab-content">{self._aiguide_section()}</div>
<div id="calltoaction-tab" class="tab-content">{self._calltoaction_section()}</div>
<div id="export-tab" class="tab-content">{self._export_section()}</div>
<div class="footer"><h3>EnteroScope Ultimate Reporter v{metadata.get('version','1.0.0')}</h3><p>University of Ghana Medical School | Brown Beckley &lt;brownbeckley94@gmail.com&gt;</p><p><a href="https://github.com/bbeckley-hub/enteroscope" target="_blank">⭐ Star us on GitHub</a></p><p>Critical Genes Tracked: Carbapenemases • ESBLs • Colistin Resistance • Tigecycline Resistance • Biofilm Formation • Efflux Pumps • Environmental Co-Selection</p><p>Generated on {metadata.get('analysis_date','Unknown')}</p></div>
</div>
<script>
    const facts = {json.dumps(self.educational_facts)};
    let factIdx = 0;
    const quoteContainer = document.getElementById('quoteContainer');
    const quoteText = document.getElementById('quoteText');
    function rotateFact() {{
        quoteContainer.style.opacity = '0';
        setTimeout(() => {{
            const fact = facts[factIdx];
            quoteText.innerHTML = '🧬 ' + fact;
            quoteContainer.style.opacity = '1';
            factIdx = (factIdx + 1) % facts.length;
        }}, 500);
    }}
    setInterval(rotateFact, 10000);
    document.addEventListener('DOMContentLoaded', rotateFact);
</script>
</body></html>"""
        return html
    
    # ---------- Section methods with enhanced research-focused descriptions ----------
    def _summary_section(self, kwargs, carb_count, bac_total):
        samples = kwargs['samples_data']
        patterns = kwargs['patterns']
        total = len(samples)
        st_unique = len(patterns.get('st_distribution', {}))
        amr_total = sum(len(db) for db in kwargs['gene_centric'].get('amr_databases', {}).values())
        vir_total = sum(len(db) for db in kwargs['gene_centric'].get('virulence_databases', {}).values())
        high_risk = len(patterns.get('high_risk_isolates', []))
        return f"""
        <div class="alert-box alert-info"><i class="fas fa-info-circle fa-2x"></i><div><h3>📊 Executive Summary – Enterobacter cloacae Complex Analysis</h3><p><strong>What this tab shows:</strong> A high-level overview of all analyzed genomes, including total sample count, MLST diversity, AMR/virulence gene frequencies, and detection of critical resistance (carbapenemases, colistin, tigecycline).</p><p><strong>Why it matters for E. cloacae research:</strong> Enterobacter cloacae complex is an emerging multidrug-resistant pathogen associated with hospital outbreaks. Understanding population structure and resistance gene frequencies is essential for surveillance, infection control, and treatment guidance.</p><p><strong>How to interpret:</strong> High numbers of unique STs suggest genetic diversity, while dominant STs indicate potential outbreak clones. Elevated carbapenemase counts (alert box) signal a high-risk CRE scenario requiring immediate intervention.</p></div></div>
        {f'<div class="alert-box alert-danger"><i class="fas fa-exclamation-triangle fa-2x"></i><div><h3>⚠️ Critical Resistance Alert</h3><p><strong>{carb_count} carbapenemase genes</strong> detected across samples. Carbapenem-resistant Enterobacteriaceae (CRE) are a critical public health threat with mortality rates >40%. Immediate infection control measures recommended.</p></div></div>' if carb_count>0 else ''}
        {f'<div class="alert-box alert-info"><i class="fas fa-globe-africa fa-2x"></i><div><h3>⚠️ Environmental Co-Selection Alert</h3><p><strong>{bac_total} environmental resistance markers</strong> detected (biocides/heavy metals). These genes can co-select for antibiotic resistance in hospital environments, promoting persistence despite disinfection.</p></div></div>' if bac_total>0 else ''}
        <h3>Key Statistics</h3>
        <div class="master-scrollable-container"><table id="summary-table" class="data-table"><thead><tr><th data-sort="string">Metric</th><th data-sort="number">Count</th><th data-sort="string">Details</th></tr></thead><tbody>
        <tr><td>Total Samples Analyzed</th><td><strong>{total}</strong></td><td>Complete genomic analysis</th></tr>
        <tr><td>Sequence Types (MLST)</th><td><strong>{st_unique}</strong></td><td>Enterobacter MLST scheme</th></tr>
        <tr><td>Total AMR Genes</th><td><strong>{amr_total}</strong></td><td>Across all AMR databases</th></tr>
        <tr><td>Carbapenemase Genes</th><td><span class="badge {'badge-critical' if carb_count>0 else 'badge-low'}">{carb_count}</span></td><td>KPC, NDM, VIM, IMP, OXA-48 types</th></tr>
        <tr><td>Virulence Genes</th><td><strong>{vir_total}</strong></td><td>Biofilm, iron uptake, toxins</th></tr>
        <tr><td>Environmental Markers</th><td><span class="badge {'badge-high' if bac_total>0 else 'badge-low'}">{bac_total}</span></td><td>Heavy metal, biocide, stress response</th></tr>
        <tr><td>High-Risk Isolates</th><td><span class="badge {'badge-critical' if high_risk>0 else 'badge-low'}">{high_risk}</span></td><td>Carbapenemase + last-resort resistance</th></tr>
        </tbody></table></div>
        """
    
    def _samples_section(self, kwargs, qc_data):
        samples = kwargs['samples_data']
        gene_centric = kwargs['gene_centric']
        sample_vir_counts = defaultdict(int)
        for db_name, genes in gene_centric.get('virulence_databases', {}).items():
            for g in genes:
                for genome in g['genomes']:
                    sample_vir_counts[genome] += 1
        html = """
        <div class="alert-box alert-info"><i class="fas fa-info-circle"></i><div><h3>🧪 Sample Overview – Per-Genome Metadata</h3><p><strong>What this tab shows:</strong> Each genome's MLST sequence type (ST), virulence gene count, and assembly quality metrics (N50, GC%, AT%).</p><p><strong>Why it matters:</strong> Linking resistance/virulence genotypes to specific isolates enables outbreak tracking, clonal comparison, and identification of hypervirulent or high-risk clones. MLST provides a portable nomenclature for global epidemiology.</p><p><strong>How to interpret:</strong> High virulence counts (e.g., >10) suggest enhanced pathogenic potential. Low N50 indicates fragmented assemblies that may miss genes. ST clustering indicates related isolates – consider if from same patient/ward/time period.</p></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('samples-tab')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-secondary" onclick="document.getElementById('search-samples').value=''; searchTable('samples-table','search-samples')"><i class="fas fa-sync"></i> Clear Search</button><button class="action-btn btn-success" onclick="exportTableToCSV('samples-table','enterobacter_samples.csv')"><i class="fas fa-download"></i> Export CSV</button></div>
        <input type="text" class="search-box" id="search-samples" onkeyup="searchTable('samples-table','search-samples')" placeholder="🔍 Search samples by any field...">
        <div class="master-scrollable-container"><table id="samples-table" class="data-table"><thead><tr><th data-sort="string">Sample</th><th data-sort="string">Sequence Type</th><th data-sort="number">Virulence Count</th><th data-sort="number">N50</th><th data-sort="number">GC%</th><th data-sort="number">AT%</th></tr></thead><tbody>
        """
        for sample, data in samples.items():
            st = data.get('mlst', {}).get('ST', 'ND')
            vcount = sample_vir_counts.get(sample, 0)
            st_disp = f"ST{st}" if st != 'ND' else 'ND'
            n50 = qc_data.get(sample, {}).get('N50', 'ND')
            if isinstance(n50, float) and n50 != 'ND':
                n50 = f"{n50:,.0f}"
            gc = qc_data.get(sample, {}).get('GC%', 'ND')
            if isinstance(gc, float):
                gc = f"{gc:.1f}"
            at = qc_data.get(sample, {}).get('AT%', 'ND')
            if isinstance(at, float):
                at = f"{at:.1f}"
            html += f'<tr><td><strong>{sample}</strong></td><td>{st_disp}</td><td>{vcount}</td><td>{n50}</td><td>{gc}%</td><td>{at}%</td></tr>'
        html += '</tbody></table></div>'
        return html
    
    def _qc_section(self, kwargs):
        qc = kwargs.get('qc_data', {})
        if not qc:
            return '<div class="alert-box alert-warning"><i class="fas fa-exclamation-circle"></i><div><h3>FASTA QC – No Data</h3><p>No FASTA QC summary file found. Assembly metrics unavailable.</p></div></div>'
        metrics = set()
        for d in qc.values():
            metrics.update(d.keys())
        metrics = sorted(metrics)
        html = '<div class="alert-box alert-info"><i class="fas fa-chart-line"></i><div><h3>🔬 FASTA Quality Control – Assembly Metrics</h3><p><strong>What this tab shows:</strong> Key assembly statistics per genome: N50 (contig length at 50% of assembly), GC%, AT%, number of contigs, total length, etc.</p><p><strong>Why it matters:</strong> Assembly quality directly affects gene detection accuracy. Poor assemblies (high contig count, low N50) may miss resistance/virulence genes or produce false positives. GC% deviation from ~55% (typical for E. cloacae) might indicate contamination.</p><p><strong>How to interpret:</strong> N50 > 50 kb indicates good assembly; N50 > 100 kb is excellent. Total length 4.5–5.5 Mbp is expected. GC% 54–56% is normal. High contig count (>500) suggests fragmented assembly – consider re-assembly or long-read sequencing.</p></div></div>'
        html += '<div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection(\'qc-tab\')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-secondary" onclick="document.getElementById(\'search-qc\').value=\'\'; searchTable(\'qc-table\',\'search-qc\')"><i class="fas fa-sync"></i> Clear Search</button><button class="action-btn btn-success" onclick="exportTableToCSV(\'qc-table\',\'fasta_qc.csv\')"><i class="fas fa-download"></i> Export CSV</button></div>'
        html += '<input type="text" class="search-box" id="search-qc" onkeyup="searchTable(\'qc-table\',\'search-qc\')" placeholder="🔍 Search sample...">'
        html += '<div class="master-scrollable-container"><table id="qc-table" class="data-table"><thead><tr><th data-sort="string" style="min-width:200px">Sample</th>'
        for m in metrics:
            html += f'<th data-sort="number">{m}</th>'
        html += '</tr></thead><tbody>'
        for sample, vals in sorted(qc.items()):
            html += f'<tr><td style="white-space:nowrap"><strong>{sample}</strong></td>'
            for m in metrics:
                v = vals.get(m, 'ND')
                if isinstance(v, float):
                    v = f"{v:,.0f}" if v > 1000 else f"{v:.2f}"
                html += f'<td>{v}</td>'
            html += '</tr>'
        html += '</tbody></table></div>'
        return html
    
    def _mlst_section(self, kwargs):
        patterns = kwargs['patterns']
        samples = kwargs['samples_data']
        st_dist = patterns.get('st_distribution', {})
        total_st = sum(st_dist.values())
        html = f"""<div class="alert-box alert-info"><i class="fas fa-code-branch"></i><div><h3>🧬 MLST – Sequence Typing and Population Structure</h3><p><strong>What this tab shows:</strong> Distribution of sequence types (STs) among isolates, based on the Enterobacter cloacae MLST scheme (7 housekeeping genes: fusA, rpoB, etc.).</p><p><strong>Why it matters:</strong> MLST reveals clonal relationships, outbreak clusters, and internationally disseminated high-risk clones (e.g., ST171, ST78, ST114, ST418). Certain STs are associated with carbapenemase acquisition and hospital spread.</p><p><strong>How to interpret:</strong> A single dominant ST suggests a common source outbreak. Multiple unique STs indicate sporadic cases or environmental diversity. Novel STs (not in PubMLST) may represent new lineages. Cross-reference STs with AMR profiles to identify high-risk clones.</p><p>Detected <strong>{len(st_dist)} unique STs</strong> among {total_st} isolates.</p></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('mlst-tab')"><i class="fas fa-print"></i> Print Section</button></div>
        <div class="master-scrollable-container"><table id="mlst-table" class="data-table"><thead><tr><th data-sort="string">ST</th><th data-sort="number">Frequency</th><th data-sort="string">Samples</th></tr></thead><tbody>"""
        for st, cnt in sorted(st_dist.items(), key=lambda x: x[1], reverse=True):
            if st == 'ND': continue
            pct = (cnt/total_st)*100 if total_st else 0
            sample_list = ', '.join([s for s,d in samples.items() if d.get('mlst',{}).get('ST')==st])
            html += f'<tr><td><strong>ST{st}</strong></td><td>{cnt} ({pct:.1f}%)</td><td>{sample_list}</td></tr>'
        html += '</tbody></table></div>'
        return html
    
    def _amr_section(self, kwargs):
        gene_centric = kwargs['gene_centric']
        all_genes = []
        for db in gene_centric.get('amr_databases', {}).values():
            all_genes.extend(db)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        filters = [('KPC','KPC'),('NDM','NDM'),('VIM','VIM'),('IMP','IMP'),('OXA-48','OXA-48'),
                   ('CTX-M','CTX-M'),('SHV','SHV'),('CMY','CMY'),('mcr','mcr'),('tet(X)','tet(X)'),
                   ('aac','aac'),('aph','aph'),('erm','erm'),('qnr','qnr'),('tet','tet'),('sul','sul')]
        html = """
        <div class="alert-box alert-info"><i class="fas fa-biohazard"></i><div><h3>💊 Antimicrobial Resistance Genes – Comprehensive AMR Profile</h3><p><strong>What this tab shows:</strong> All detected AMR genes across multiple databases (AMRfinder, CARD, ResFinder, etc.), including carbapenemases, ESBLs, colistin resistance (mcr), and tigecycline resistance (tet(X), efflux pumps).</p><p><strong>Why it matters for E. cloacae:</strong> This species is a major CRE pathogen. Detecting carbapenemases (KPC, NDM, OXA-48) guides therapy (e.g., ceftazidime-avibactam vs. aztreonam-avibactam). Colistin resistance (mcr) or tigecycline resistance predicts treatment failure of last-resort drugs.</p><p><strong>How to interpret:</strong> Sort by frequency to see common resistance determinants. Click on filter buttons to focus on high-priority genes. The "Genomes" column lists which isolates carry the gene – use the second search box to highlight specific genomes across the table.</p></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('amr-tab')"><i class="fas fa-print"></i> Print Section</button></div>
        <input type="text" class="search-box" id="search-amr-gene" onkeyup="searchTable('amr-table','search-amr-gene')" placeholder="🔍 Search AMR genes...">
        <input type="text" class="search-box" id="search-amr-genome" onkeyup="highlightGenome('amr-table','search-amr-genome')" placeholder="🔍 Highlight genomes containing the text">
        <div class="action-buttons">"""
        for display, term in filters:
            html += f'<button class="action-btn btn-warning" onclick="document.getElementById(\'search-amr-gene\').value=\'{term}\'; searchTable(\'amr-table\',\'search-amr-gene\')">{display}</button>'
        html += '<button class="action-btn btn-primary" onclick="document.getElementById(\'search-amr-gene\').value=\'\'; searchTable(\'amr-table\',\'search-amr-gene\')">Clear</button>'
        html += '</div><div class="master-scrollable-container"><table id="amr-table" class="data-table"><thead><tr><th data-sort="string">Gene</th><th data-sort="string">Database</th><th data-sort="number">Frequency</th><th class="col-genomes" data-sort="string">Genomes</th></tr></thead><tbody>'
        for g in all_genes:
            tags = ''.join(f'<span class="genome-tag">{gen}</span>' for gen in g['genomes'])
            disp = f"<strong>{g['gene']}</strong>" + (' 🔥' if g['category']=='Carbapenemases' else '')
            html += f'<tr><td>{disp}</td><td>{g["database"]}</td><td>{g["frequency_display"]}</td><td class="col-genomes"><div class="genome-list">{tags}</div></td></tr>'
        html += '</tbody></table></div>'
        return html
    
    def _virulence_section(self, kwargs):
        gene_centric = kwargs['gene_centric']
        all_genes = []
        for db in gene_centric.get('virulence_databases', {}).values():
            all_genes.extend(db)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        filters = ['ompA', 'csu', 'bfmR', 'bfmS', 'bap', 'pga', 'csg', 'fim', 'tss', 'hcp', 'vgrG', 'hly', 'ent', 'iuc', 'iro']
        html = """
        <div class="alert-box alert-info"><i class="fas fa-virus"></i><div><h3>🦠 Virulence Factors – Pathogenicity Determinants</h3><p><strong>What this tab shows:</strong> Virulence-associated genes including biofilm formation (ompA, pgaABCD, csg), adhesion, type VI secretion systems (T6SS), iron acquisition (enterobactin, aerobactin), and toxins.</p><p><strong>Why it matters:</strong> Biofilm genes enable persistent infections on medical devices (catheters, ventilators). Iron acquisition systems enhance survival in host tissues. T6SS mediates interbacterial competition and may modulate host immune response.</p><p><strong>How to interpret:</strong> High prevalence of biofilm genes suggests these isolates are capable of forming recalcitrant biofilms – important for infection control. Search/click filters to check specific virulence mechanisms. Cross-reference with AMR genes to identify isolates with both high resistance and high virulence potential.</p></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('virulence-tab')"><i class="fas fa-print"></i> Print Section</button></div>
        <input type="text" class="search-box" id="search-vir-gene" onkeyup="searchTable('vir-table','search-vir-gene')" placeholder="🔍 Search virulence genes...">
        <input type="text" class="search-box" id="search-vir-genome" onkeyup="highlightGenome('vir-table','search-vir-genome')" placeholder="🔍 Highlight genomes">
        <div class="action-buttons">"""
        for f in filters:
            html += f'<button class="action-btn btn-success" onclick="document.getElementById(\'search-vir-gene\').value=\'{f}\'; searchTable(\'vir-table\',\'search-vir-gene\')">{f}</button>'
        html += '<button class="action-btn btn-primary" onclick="document.getElementById(\'search-vir-gene\').value=\'\'; searchTable(\'vir-table\',\'search-vir-gene\')">Clear</button>'
        html += '</div><div class="master-scrollable-container"><table id="vir-table" class="data-table"><thead><tr><th data-sort="string">Gene</th><th data-sort="string">Database</th><th data-sort="number">Frequency</th><th class="col-genomes" data-sort="string">Genomes</th></tr></thead><tbody>'
        for g in all_genes:
            tags = ''.join(f'<span class="genome-tag">{gen}</span>' for gen in g['genomes'])
            html += f'<tr><td><strong>{g["gene"]}</strong></td><td>{g["database"]}</td><td>{g["frequency_display"]}</td><td class="col-genomes"><div class="genome-list">{tags}</div></td></tr>'
        html += '</tbody></table></div>'
        return html
    
    def _bacmet_section(self, kwargs):
        gene_centric = kwargs['gene_centric']
        all_genes = []
        for db in gene_centric.get('bacmet_databases', {}).values():
            all_genes.extend(db)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        filters = ['qac', 'cep', 'form', 'mer', 'ars', 'cop', 'sil', 'chr', 'cad', 'znt', 'czc', 'nik', 'ade', 'mex']
        html = """
        <div class="alert-box alert-info"><i class="fas fa-flask"></i><div><h3>🧪 Biocide & Heavy Metal Resistance – Environmental Co-Selection Markers</h3><p><strong>What this tab shows:</strong> Genes conferring resistance to disinfectants (qac, cepA, form), heavy metals (mercury – mer, arsenic – ars, copper – cop, silver – sil), and other environmental stressors.</p><p><strong>Why it matters:</strong> These markers often co-locate with AMR genes on mobile genetic elements (plasmids, integrons). Their presence indicates potential co-selection of antibiotic resistance in hospital environments exposed to biocides and heavy metals (e.g., from disinfectants, plumbing, medical devices).</p><p><strong>How to interpret:</strong> High prevalence of qac genes suggests reduced susceptibility to quaternary ammonium compounds (common disinfectants). Mercury/arsenic resistance may indicate historical or ongoing heavy metal exposure. Search for specific markers to assess disinfection efficacy risk.</p></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('bacmet-tab')"><i class="fas fa-print"></i> Print Section</button></div>
        <input type="text" class="search-box" id="search-bac-gene" onkeyup="searchTable('bac-table','search-bac-gene')" placeholder="🔍 Search BACMET genes...">
        <input type="text" class="search-box" id="search-bac-genome" onkeyup="highlightGenome('bac-table','search-bac-genome')" placeholder="🔍 Highlight genomes">
        <div class="action-buttons">"""
        for f in filters:
            html += f'<button class="action-btn btn-warning" onclick="document.getElementById(\'search-bac-gene\').value=\'{f}\'; searchTable(\'bac-table\',\'search-bac-gene\')">{f}</button>'
        html += '<button class="action-btn btn-primary" onclick="document.getElementById(\'search-bac-gene\').value=\'\'; searchTable(\'bac-table\',\'search-bac-gene\')">Clear</button>'
        html += '</div><div class="master-scrollable-container"><table id="bac-table" class="data-table"><thead><tr><th data-sort="string">Gene</th><th data-sort="string">Database</th><th data-sort="number">Frequency</th><th class="col-genomes" data-sort="string">Genomes</th></tr></thead><tbody>'
        for g in all_genes:
            tags = ''.join(f'<span class="genome-tag">{gen}</span>' for gen in g['genomes'])
            html += f'<tr><td><strong>{g["gene"]}</strong></td><td>{g["database"]}</td><td>{g["frequency_display"]}</td><td class="col-genomes"><div class="genome-list">{tags}</div></td></tr>'
        html += '</tbody></table></div>'
        return html
    
    def _plasmid_section(self, kwargs):
        plasmid = kwargs.get('plasmid_analysis', {})
        genes = plasmid.get('plasmid_frequencies', [])
        html = """<div class="alert-box alert-info"><i class="fas fa-dna"></i><div><h3>🧬 Plasmid Replicons – Mobile Genetic Elements</h3><p><strong>What this tab shows:</strong> Plasmid replicon markers (e.g., IncL/M, IncFII, IncFIB, IncHI2, Col-like) and their distribution among isolates.</p><p><strong>Why it matters for E. cloacae:</strong> Epidemic carbapenemase genes (blaOXA-48, blaKPC, blaNDM) are often carried on specific plasmids: IncL/M frequently carries blaOXA-48, IncFII/IncFIB often carry blaKPC/blaCTX-M, and IncHI2 may carry mcr or NDM. Identifying plasmid replicons helps predict resistance gene mobility and outbreak potential.</p><p><strong>How to interpret:</strong> High-frequency replicons indicate successful plasmid backbones spreading in the population. Presence of IncL/M + blaOXA-48 suggests a transferable carbapenem resistance threat. Multiple replicons in the same isolate may indicate plasmid accumulation. If no data, run PlasmidFinder separately.</p></div></div>"""
        if not genes:
            return html + '<p><strong>📭 No PlasmidFinder data found.</strong> This tab requires a PlasmidFinder summary report. Please run EnteroScope with plasmid screening enabled.</p>'
        html += '<div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection(\'plasmid-tab\')"><i class="fas fa-print"></i> Print Section</button></div>'
        html += '<input type="text" class="search-box" id="search-plasmid-gene" onkeyup="searchTable(\'plasmid-table\',\'search-plasmid-gene\')" placeholder="🔍 Search plasmids...">'
        html += '<input type="text" class="search-box" id="search-plasmid-genome" onkeyup="highlightGenome(\'plasmid-table\',\'search-plasmid-genome\')" placeholder="🔍 Highlight genomes">'
        html += '<div class="master-scrollable-container"><table id="plasmid-table" class="data-table"><thead><tr><th data-sort="string">Marker</th><th data-sort="string">Category</th><th data-sort="number">Frequency</th><th class="col-genomes" data-sort="string">Genomes</th></tr></thead><tbody>'
        for g in genes:
            tags = ''.join(f'<span class="genome-tag">{gen}</span>' for gen in g['genomes'])
            html += f'<tr><td><strong>{g["plasmid_marker"]}</strong></td><td>{g["category"]}</td><td>{g["frequency_display"]}</td><td class="col-genomes"><div class="genome-list">{tags}</div></td></tr>'
        html += '</tbody></table></div>'
        return html
    
    def _patterns_section(self, kwargs):
        patterns = kwargs['patterns']
        high_risk = patterns.get('high_risk_isolates', [])
        gene_centric = kwargs.get('gene_centric', {})
        samples_data = kwargs.get('samples_data', {})
        total_genomes = len(samples_data)
        
        # Build a per‑genome set of all detected genes (all categories)
        sample_genes = {}
        for db_type in ['amr_databases', 'virulence_databases', 'bacmet_databases', 'plasmid_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for g in genes:
                    for genome in g['genomes']:
                        sample_genes.setdefault(genome, set()).add(g['gene'])
        
        # Compute co‑occurrence counts (A,B) with A < B
        cooc = defaultdict(int)
        genes_list = list(set().union(*sample_genes.values()))
        for g in genes_list:
            for h in genes_list:
                if g < h:
                    count = 0
                    for genome, geneset in sample_genes.items():
                        if g in geneset and h in geneset:
                            count += 1
                    if count > 0:
                        cooc[(g, h)] = count
        
        # Get top 500 by count
        top_cooc = sorted(cooc.items(), key=lambda x: x[1], reverse=True)[:500]
        
        # Generate HTML for co‑occurrence table
        cooc_html = ""
        if top_cooc:
            cooc_html = """
            <hr style="margin:40px 0 20px 0;">
            <div class="alert-box alert-info">
                <i class="fas fa-link fa-2x"></i>
                <div>
                    <h3>🔗 Gene Co‑occurrence (Top 500 Pairs)</h3>
                    <p><strong>What this shows:</strong> Pairs of genes that frequently appear together in the same genome. High co‑occurrence may indicate physical linkage (same plasmid, integron, or genomic island) or shared selection pressure.</p>
                    <p><strong>How to use:</strong> Use the search box below to filter by gene name. Sort by count to find the most common associations (e.g., carbapenemase + efflux pump, or mcr + Inc plasmid).</p>
                </div>
            </div>
            <div class="action-buttons">
                <button class="action-btn btn-primary" onclick="printSection('patterns-tab')"><i class="fas fa-print"></i> Print Section</button>
                <button class="action-btn btn-success" onclick="exportTableToCSV('cooc-table','gene_cooccurrence.csv')"><i class="fas fa-download"></i> Export CSV</button>
            </div>
            <input type="text" class="search-box" id="search-cooc" onkeyup="searchTable('cooc-table','search-cooc')" placeholder="🔍 Search genes in co‑occurrence table...">
            <div class="master-scrollable-container">
                <table id="cooc-table" class="data-table">
                    <thead>
                        <tr>
                            <th data-sort="string">Gene A</th>
                            <th data-sort="string">Gene B</th>
                            <th data-sort="number">Co‑occurrence Count</th>
                            <th data-sort="number">Percentage (%)</th>
                        </tr>
                    </thead>
                    <tbody>
            """
            for (gA, gB), cnt in top_cooc:
                pct = (cnt / total_genomes) * 100 if total_genomes else 0
                cooc_html += f"<tr><td><strong>{gA}</strong></td><td><strong>{gB}</strong></td><td>{cnt}</td><td>{pct:.1f}%</td></tr>"
            cooc_html += """
                    </tbody>
                </table>
            </div>
            """
        else:
            cooc_html = "<p>No co‑occurrence data available (fewer than 2 genes detected).</p>"
        
        # High‑risk isolates table (fixed categorization ensures fimI is not shown as carbapenemase)
        risk_html = """<div class="alert-box alert-info"><i class="fas fa-project-diagram"></i><div><h3>⚠️ High‑Risk Isolates – Carbapenemase + Last‑Line Resistance</h3><p><strong>What this tab shows:</strong> Isolates that co‑carry carbapenemase genes AND resistance to at least one last‑resort antibiotic (colistin, tigecycline).</p><p><strong>Why it matters:</strong> These "pan-resistant" isolates leave virtually no treatment options. Their emergence signals a critical threat requiring urgent public health response, infection control, and alternative therapies (e.g., cefiderocol, phage therapy).</p><p><strong>How to interpret:</strong> Each row lists the sample, ST, and specific carbapenemase/colistin/tigecycline genes. Cross-reference with MLST tab to see if high‑risk isolates share the same ST (possible outbreak). Immediate action: isolate patients, enhanced environmental cleaning, and contact precautions.</p></div></div>"""
        if high_risk:
            risk_html += '<div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection(\'patterns-tab\')"><i class="fas fa-print"></i> Print Section</button></div>'
            risk_html += '<input type="text" class="search-box" id="search-patterns" onkeyup="searchTable(\'crab-table\',\'search-patterns\')" placeholder="🔍 Search high‑risk isolates...">'
            risk_html += '<div class="master-scrollable-container"><table id="crab-table" class="data-table"><thead><tr><th data-sort="string">Sample</th><th data-sort="string">ST</th><th data-sort="string">Carbapenemases</th><th data-sort="string">Colistin Resistance</th><th data-sort="string">Tigecycline Resistance</th></tr></thead><tbody>'
            for c in high_risk:
                risk_html += f'<tr><td><strong>{c["sample"]}</strong></td><td>{c.get("st","ND")}</td><td>{",".join(c["carbapenemases"])}</td><td>{",".join(c["colistin_resistance"])}</td><td>{",".join(c["tigecycline_resistance"])}</td></tr>'
            risk_html += '</tbody></table></div>'
        else:
            risk_html += '<p>✅ No high‑risk isolates (carbapenemase + colistin/tigecycline resistance) detected in this dataset.</p>'
        
        return risk_html + cooc_html
    
    def _databases_section(self, kwargs):
        stats = kwargs['gene_centric'].get('database_stats', {})
        html = """<div class="alert-box alert-info"><i class="fas fa-database"></i><div><h3>📚 Database Metrics – Gene Coverage Summary</h3><p><strong>What this tab shows:</strong> Per-database statistics: number of unique genes, total occurrences (sum of gene counts across genomes), and critical genes (carbapenemases, ESBLs, colistin, tigecycline).</p><p><strong>Why it matters:</strong> Different databases have varying coverage and specificity. AMRfinder provides FDA-ARGOS curated genes, CARD includes resistance ontology, ResFinder focuses on acquired resistance, and BacMet covers biocides/metals.</p><p><strong>How to interpret:</strong> Larger unique gene counts indicate broader coverage. High total occurrences suggest widespread resistance determinants. Compare databases to identify genes found by multiple tools (higher confidence) vs. unique hits (may need validation).</p></div></div>"""
        html += '<div class="master-scrollable-container"><table class="data-table"><thead><tr><th data-sort="string">Database</th><th data-sort="number">Unique Genes</th><th data-sort="number">Total Occurrences</th><th data-sort="number">Critical Genes</th></tr></thead><tbody>'
        for db, d in stats.items():
            html += f'<tr><td><strong>{db.upper()}</strong></td><td>{d["total_genes"]}</td><td>{d["total_occurrences"]}</td><td>{d["critical_genes"]}</td></tr>'
        html += '</tbody></table></div>'
        return html
    
    def _aiguide_section(self):
        return """
        <div class="alert-box alert-info"><i class="fas fa-robot fa-2x"></i><div><h3>🤖 AI Assistant Integration Guide</h3><p><strong>What this tab shows:</strong> Instructions for using large language models (ChatGPT, Claude, Gemini, DeepSeek) to analyze your Enterobacter cloacae genome data.</p><p><strong>Why it matters:</strong> AI can help interpret complex resistance patterns, suggest treatment correlations, generate hypotheses, and draft reports – accelerating research and clinical decision-making.</p><p><strong>How to use:</strong> Upload this HTML file or the companion JSON file to your preferred AI assistant. Then ask questions such as:</p><ul><li>Which samples carry carbapenemase genes?</li><li>List all isolates with high biofilm‑related gene counts.</li><li>What are the most common biocide resistance genes in this dataset?</li><li>Show me ST-plasmid associations.</li><li>Which high-risk isolates require immediate infection control?</li><li>Summarize the AMR profile of sample X.</li></ul><p>The AI can also generate custom tables, maps, and even draft manuscripts based on these data.</p></div></div>
        """
    
    def _calltoaction_section(self):
        return """
        <div class="alert-box alert-info"><i class="fas fa-globe fa-2x"></i><div><h3>🌍 Call to Action – Combating AMR Together</h3><p><strong>Global burden:</strong> Antimicrobial resistance (AMR) kills an estimated 1.27 million people annually, with carbapenem-resistant Enterobacterales (CRE) classified as a critical priority pathogen by WHO. <em>Enterobacter cloacae</em> complex is a major contributor to hospital-acquired CRE infections worldwide.</p><p><strong>Concrete actions you can take:</strong></p><ul><li><i class="fas fa-share-alt"></i> <strong>Share data to public databases:</strong> Submit genomes and resistance profiles to NCBI (BioProject), PubMLST, and PLSDB to enable global surveillance.</li><li><i class="fas fa-chart-line"></i> <strong>Implement local surveillance:</strong> Establish routine genomic surveillance of Enterobacter cloacae isolates in your hospital or region to detect emerging resistance early.</li><li><i class="fas fa-hand-holding-heart"></i> <strong>Adopt infection control measures:</strong> Use contact precautions, enhanced environmental disinfection, and antimicrobial stewardship to limit spread of CRE.</li><li><i class="fab fa-github"></i> <strong>Contribute to open science:</strong> Improve EnteroScope, report bugs, suggest features, or fork the repository on GitHub. Your contributions benefit the entire research community.</li><li><i class="fas fa-envelope"></i> <strong>Collaborate with us:</strong> Reach out for collaborations, data sharing, or joint grant applications. We welcome partnerships with clinicians, epidemiologists, and bioinformaticians.</li></ul><div style="text-align:center; margin:30px 0;"><i class="fas fa-star" style="font-size:3em; color:#ffc107;"></i><br/><a href="https://github.com/bbeckley-hub/enteroscope" target="_blank" style="font-size:1.2em;">⭐ Star us on GitHub ⭐</a><br/><br/><p><strong>📧 Email:</strong> <a href="mailto:brownbeckley94@gmail.com">brownbeckley94@gmail.com</a><br/><strong>🔗 GitHub:</strong> <a href="https://github.com/bbeckley-hub/enteroscope" target="_blank">github.com/bbeckley-hub/enteroscope</a></p></div></div></div>
        """
    
    def _export_section(self):
        return """
        <div class="alert-box alert-info"><i class="fas fa-download"></i><div><h3>📥 Export Data – CSV and JSON Formats</h3><p><strong>What this tab shows:</strong> Buttons to download tables (CSV) for further analysis in Excel, R, or Python, and the full JSON report for programmatic access or AI upload.</p><p><strong>Why it matters:</strong> Enables custom plots, statistical analyses, metadata integration, and sharing with collaborators who prefer tabular data.</p><p><strong>How to use:</strong> Click any CSV button to download the current filtered/sorted table. The JSON file contains the complete integrated dataset (use the link after report generation).</p></div></div>
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('samples-table','enterobacter_samples.csv')">Samples CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('amr-table','amr_genes.csv')">AMR CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('vir-table','virulence_genes.csv')">Virulence CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('bac-table','bacmet_genes.csv')">Bacmet CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('plasmid-table','plasmid_markers.csv')">Plasmid CSV</button>
            <button class="action-btn btn-success" onclick="location.href='enteroscope_ultimate_report.json'">Download JSON (full data)</button>
        </div>
        """

# =============================================================================
# MAIN REPORTER CLASS 
# =============================================================================
class EnteroUltimateReporter:
    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "ENTERO_ULTIMATE_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = EnteroUltimateHTMLParser()
        self.analyzer = EnteroDataAnalyzer()
        self.html_generator = EnteroHTMLGenerator(self.analyzer)
        self.metadata = {
            "tool_name": "EnteroScope Ultimate Reporter",
            "version": "1.0.0",
            "author": "Brown Beckley <brownbeckley94@gmail.com>",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir)
        }
    
    def find_html_files(self) -> Dict[str, List[Path]]:
        print("🔍 Searching for EnteroScope HTML reports...")
        html_files = {
            'mlst': [], 'amrfinder': [], 'abricate': defaultdict(list),
            'plasmidfinder': [], 'qc': []
        }
        all_html = list(self.input_dir.glob("**/*.html"))
        print(f"  📁 Found {len(all_html)} HTML files")
        for f in all_html:
            name = f.name.lower()
            if 'mlst_summary' in name:
                html_files['mlst'].append(f)
            elif 'amrfinder' in name and 'summary' in name:
                html_files['amrfinder'].append(f)
            elif 'plasmidfinder' in name and 'summary' in name:
                html_files['plasmidfinder'].append(f)
            elif 'fasta_qc_summary' in name:
                html_files['qc'].append(f)
            else:
                matched = False
                for db_key in self.parser.db_name_mapping.keys():
                    if db_key in name and 'plasmidfinder' not in db_key and 'amrfinder' not in db_key:
                        html_files['abricate'][db_key].append(f)
                        matched = True
                        break
                if not matched and any(db in name for db in ['card','resfinder','argannot','megares','bacmet2','vfdb','ecoli_vf','ecoh','ncbi']):
                    for db in ['card','resfinder','argannot','megares','bacmet2','vfdb','ecoli_vf','ecoh','ncbi']:
                        if db in name:
                            html_files['abricate'][f'enteroscope_{db}'].append(f)
                            break
        print(f"  ✅ MLST: {len(html_files['mlst'])} | AMRfinder: {len(html_files['amrfinder'])} | PlasmidFinder: {len(html_files['plasmidfinder'])} | QC: {len(html_files['qc'])} | ABRicate DBs: {len(html_files['abricate'])}")
        return html_files
    
    def integrate_all_data(self, html_files: Dict[str, List[Path]]) -> Dict[str, Any]:
        print("\n🔗 Integrating data...")
        integrated = {'metadata': self.metadata, 'samples': {}, 'patterns': {}, 'gene_centric': {}, 'plasmid_analysis': {}, 'qc_data': {}}
        if html_files['qc']:
            integrated['qc_data'] = self.parser.parse_qc_report(html_files['qc'][0])
        mlst_data = self.parser.parse_mlst_report(html_files['mlst'][0]) if html_files['mlst'] else {}
        total_samples = len(mlst_data)
        amr_by_sample, amr_gene_freq = {}, {}
        if html_files['amrfinder']:
            amr_by_sample, amr_gene_freq = self.parser.parse_amrfinder_report(html_files['amrfinder'][0], total_samples)
        abricate_by_sample = defaultdict(dict)
        abricate_gene_freq = {}
        for db_key, files in html_files['abricate'].items():
            if files:
                db_name = self.parser.db_name_mapping.get(db_key, db_key)
                by_sample, freq = self.parser.parse_abricate_database_report(files[0], total_samples)
                for s, genes in by_sample.items():
                    abricate_by_sample[s][db_name] = genes
                abricate_gene_freq[db_name] = freq
        plasmid_by_sample, plasmid_gene_freq = {}, {}
        if html_files['plasmidfinder']:
            plasmid_by_sample, plasmid_gene_freq = self.parser.parse_plasmidfinder_report(html_files['plasmidfinder'][0], total_samples)
        all_samples = set(mlst_data.keys()) | set(amr_by_sample.keys()) | set(abricate_by_sample.keys()) | set(plasmid_by_sample.keys()) | set(integrated['qc_data'].keys())
        total_samples = len(all_samples)
        print(f"📊 Found {total_samples} unique samples")
        for sample in all_samples:
            virulence = []
            for db in ['vfdb', 'ecoli_vf']:
                if db in abricate_by_sample.get(sample, {}):
                    virulence.extend(abricate_by_sample[sample][db])
            integrated['samples'][sample] = {
                'mlst': mlst_data.get(sample, {'ST':'ND', 'Allele_Profile':'ND'}),
                'amr_genes': amr_by_sample.get(sample, []),
                'virulence_genes': list(set(virulence)),
                'environmental_genes': abricate_by_sample.get(sample, {}).get('bacmet2', []),
                'plasmid_genes': plasmid_by_sample.get(sample, [])
            }
        integrated['gene_frequencies'] = {'amrfinder': amr_gene_freq, 'abricate': abricate_gene_freq}
        if plasmid_gene_freq:
            integrated['gene_frequencies']['plasmidfinder'] = plasmid_gene_freq
        integrated['gene_centric'] = self.analyzer.create_gene_centric_tables(integrated, total_samples)
        integrated['patterns'] = self.analyzer.create_cross_genome_patterns(integrated, total_samples)
        integrated['plasmid_analysis'] = self.analyzer.create_plasmid_analysis(integrated, total_samples)
        return integrated
    
    def generate_json_report(self, data: Dict[str, Any]) -> Path:
        print("\n📝 Generating JSON report...")
        output_file = self.output_dir / "enteroscope_ultimate_report.json"
        def convert(obj):
            if isinstance(obj, dict):
                return {k: convert(v) for k, v in obj.items()}
            elif isinstance(obj, (list, tuple, set)):
                return [convert(i) for i in obj]
            elif isinstance(obj, defaultdict):
                return convert(dict(obj))
            else:
                return obj
        serializable = convert(data)
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(serializable, f, indent=2, default=str)
        print(f"    ✅ JSON saved: {output_file}")
        return output_file
    
    def generate_csv_reports(self, data: Dict[str, Any]):
        print("\n📊 Generating CSV reports...")
        qc = data.get('qc_data', {})
        rows = []
        for s, d in data['samples'].items():
            st = d['mlst']['ST']
            vcount = len(d['virulence_genes'])
            n50 = qc.get(s, {}).get('N50', 'ND')
            gc = qc.get(s, {}).get('GC%', 'ND')
            at = qc.get(s, {}).get('AT%', 'ND')
            rows.append({'Sample': s, 'ST': st, 'Virulence_Count': vcount, 'N50': n50, 'GC%': gc, 'AT%': at})
        pd.DataFrame(rows).to_csv(self.output_dir / "enterobacter_samples.csv", index=False)
        amr_rows = []
        for db, genes in data['gene_centric'].get('amr_databases', {}).items():
            for g in genes:
                amr_rows.append({'Gene':g['gene'], 'Database':g['database'], 'Frequency':g['frequency_display'], 'Genomes':';'.join(g['genomes'])})
        if amr_rows:
            pd.DataFrame(amr_rows).to_csv(self.output_dir / "amr_genes.csv", index=False)
        vir_rows = []
        for db, genes in data['gene_centric'].get('virulence_databases', {}).items():
            for g in genes:
                vir_rows.append({'Gene':g['gene'], 'Database':g['database'], 'Frequency':g['frequency_display'], 'Genomes':';'.join(g['genomes'])})
        if vir_rows:
            pd.DataFrame(vir_rows).to_csv(self.output_dir / "virulence_genes.csv", index=False)
        bac_rows = []
        for db, genes in data['gene_centric'].get('bacmet_databases', {}).items():
            for g in genes:
                bac_rows.append({'Gene':g['gene'], 'Frequency':g['frequency_display'], 'Genomes':';'.join(g['genomes'])})
        if bac_rows:
            pd.DataFrame(bac_rows).to_csv(self.output_dir / "bacmet_genes.csv", index=False)
        plas_rows = []
        for db, genes in data.get('plasmid_analysis', {}).get('plasmid_databases', {}).items():
            for g in genes:
                plas_rows.append({'Marker':g['plasmid_marker'], 'Frequency':g['frequency_display'], 'Genomes':';'.join(g['genomes'])})
        if plas_rows:
            pd.DataFrame(plas_rows).to_csv(self.output_dir / "plasmid_markers.csv", index=False)
    
    def run(self):
        print("="*80)
        print("🧠 EnteroScope Ultimate Reporter v1.0.0 (Enterobacter cloacae Complex)")
        print("="*80)
        html_files = self.find_html_files()
        if not any(html_files.values()):
            print("❌ No HTML files found!")
            return False
        data = self.integrate_all_data(html_files)
        if not data:
            return False
        self.generate_json_report(data)
        self.generate_csv_reports(data)
        self.html_generator.generate_main_report(data, self.output_dir)
        print("\n✅ Analysis complete! Open enteroscope_ultimate_report.html in your browser.")
        return True

def main():
    parser = argparse.ArgumentParser(description='EnteroScope Ultimate Reporter for Enterobacter cloacae Complex')
    parser.add_argument('-i', '--input-dir', required=True, help='Directory with EnteroScope HTML summary reports')
    args = parser.parse_args()
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        print(f"❌ Input directory not found: {input_dir}")
        sys.exit(1)
    reporter = EnteroUltimateReporter(input_dir)
    success = reporter.run()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()