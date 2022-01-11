#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <set>
#include <math.h>      
#include <random>

using namespace std;

class Contig
{
public:
	string chr;
	string moltype = "";
	vector<int> seq;
	vector<int> repeat;
	string id;
};

class SGD
{
public:
	int chr = -1;
	int start = -1;
	int end = -1;
	int strand = -1;
	int cds_start = -1;
	int cds_end = -1;
	bool splice = false;
	string annotation_type;
	string ID;
	string Name;
	string Parent;
	string orf_classification;
	int tifseqs_hits0 = 0;
	int tifseqs_hits1 = 0;
	bool operator<(const SGD& ann)
	{
		if (ann.chr == this->chr)
		{
			return this->start < ann.start;
		}
		return this->chr < ann.chr;
	}
};


class OverlapHit
{
public:
	string annotation;
	string classification;
	int start_pos = -1;
	int end_pos = -1;
};

class SharedTranscript
{
public:
	string orf_id;
	string relation;
	string shared_sgd;
	string shared_sgd_class;
};

class Tifseq
{
public:
	int chr = -1;
	int start = -1;
	int end = -1;
	int strand = -1;
	int hits0 = -1;
	int hits1 = -1;
	vector<int> covered_orfs;
	vector<int> covered_genes;
	bool operator<(const Tifseq& ann)
	{
		if (ann.chr == this->chr)
		{
			return this->start < ann.start;
		}
		return this->chr < ann.chr;
	}
};

class Orf
{
public:
	string annotation;
	int id;
	string orf_id;
	string orf_class="None";
	bool annotated = false;
	int strand = -1;
	int start_pos = -1;
	int end_pos = -1;
	int contig = -1;
	int splice = 0;
	vector<int> seq;
	string nucseq;
	string is_gene = "X";
	string splice_gene = "X";
	string up_extend_gene = "X";
	vector<OverlapHit> overlaps;
	string overlaps_gene = "X";
	int overlap_gene_relative_frame = -1;
	int overlap_gene_length = 0;
	int tifseqs_hits0 = 0;
	int tifseqs_hits1 = 0;
	string transcript_relation;
	vector<string> coupled_genes;
	vector<int> frame_hits;
	vector<int> sum_hits;
	vector<int> riboseq_reads;
	vector<int> scrambled_reads;
	int sum_mapped_reads = 0;
	vector<double> align_evalues;
	vector<vector<int>> align_seq;
	bool homolog = false;
	vector<int> frame_seq;
	vector<int> align_index;
	int in_start = -1;
	int in_end = -1;
	vector<int> alt_start_pos;
	vector<int> alt_end_pos;
	vector<int> alt_contig;
	int ref_start = -1;
	int ref_end = -1;
	int ref_contig = -1;
	vector<int> good_align;
	vector<vector<int>> focal_align_seq;
	vector<double> align_identity;
	vector<int> overlap;
	vector<vector<int>> match_frame;
	vector<int> homo_length;
	vector<vector<int>> diffs;
	vector<int> syn_diffs;
	vector<int> nonsyn_diffs;
	vector<double> syn_diffs_exp;
	vector<double> nonsyn_diffs_exp;
	vector<vector<int>> alt_seqs;
	vector<string> alt_aaseq;
	vector<vector<int>> focal_flank_align;
	vector<vector<int>> flank_align;
	vector<int> blast_align_length;
	vector<double> blast_evalue;
	vector<int> use_blast_align;
	vector<int> homology_validated;
	int block_id = -1;
	int within_block_id = -1;
	int syn_vars = 0;
	int nonsyn_vars = 0;
	double syn_exp = 0;
	double nonsyn_exp = 0;
	double nuc_diverse = 0;
	vector<string> focal_align_aaseq;
	vector<string> alt_align_aaseq;
	int repeat_count = 0;
	double best_homo_blast =1000;

	/*int overlaps_gene_id = -1;
	string is_dub = "X";
	vector<int> exon_starts;
	vector<int> exon_end;
	int splice = 0;
	int focal_end = -1;
	int focal_contig = -1;
	*/
	string aaseq;
	double coding_score = 0;
	
	/*
	int mapped_peptides = 0;
	int potential_starts = 0;
	int earliest_cys_tyr = 0;
	vector<double> ms_qval;*/
	int has_atg = 0;
	int absent_count = 0;
	int no_atg = 0;
	int has_stop = 0;
	int inter_stop = 0;
	int conserved_orf = 0;
	/*vector<double> rpbs_bayes;
	vector<double> price_pval;
	vector<double> riborf_pval;
	bool de_novo = false;
	vector<int> outgroup;
	vector<vector<int>> mapped_reads;
	vector<int> preseq;
	vector<int> base_nucs;
	vector<int> sub_nucs;
	int base_nucs_count = 0;
	int sub_nucs_count = 0;
	vector<int> start_conserved;
	int up_gene = -1;
	int down_gene = -1;*/
	bool start_codon_overlap = false;
	bool stop_codon_overlap = false;
	
	string up_gene_sys = "";
	string down_gene_sys = "";
	int up_gene_dist = 0;
	int down_gene_dist = 0;

	Orf(const int SPECIES_COUNT)
	{
		frame_hits = vector<int>(3);
		sum_hits = vector<int>(3);
		align_evalues = vector<double>(SPECIES_COUNT, 100);

		align_seq = vector<vector<int>>(SPECIES_COUNT);
		alt_start_pos = vector<int>(SPECIES_COUNT);
		alt_end_pos = vector<int>(SPECIES_COUNT);
		alt_contig = vector<int>(SPECIES_COUNT);

		focal_align_seq = vector<vector<int>>(SPECIES_COUNT);
		flank_align = vector<vector<int>>(SPECIES_COUNT);
		focal_flank_align = vector<vector<int>>(SPECIES_COUNT);
		use_blast_align = vector<int>(SPECIES_COUNT);
		homology_validated = vector<int>(SPECIES_COUNT);

		good_align = vector<int>(SPECIES_COUNT);
		overlap = vector<int>(SPECIES_COUNT);
		match_frame = vector<vector<int>>(SPECIES_COUNT);
		diffs = vector<vector<int>>(SPECIES_COUNT, vector<int>(3));
		blast_align_length = vector<int>(SPECIES_COUNT);
		blast_evalue = vector<double>(SPECIES_COUNT);
		homo_length = vector<int>(SPECIES_COUNT);
		syn_diffs = vector<int>(SPECIES_COUNT);
		nonsyn_diffs = vector<int>(SPECIES_COUNT);
		syn_diffs_exp = vector<double>(SPECIES_COUNT);
		nonsyn_diffs_exp = vector<double>(SPECIES_COUNT);
		align_identity = vector<double>(SPECIES_COUNT);
		//outgroup = vector<int>(2, -1);
		alt_aaseq = vector<string>(SPECIES_COUNT);
		alt_seqs = vector<vector<int>>(SPECIES_COUNT);
		focal_align_aaseq = vector<string>(SPECIES_COUNT);
		alt_align_aaseq = vector<string>(SPECIES_COUNT);
	}
	bool operator<(const Orf& ann)
	{
		if (ann.contig == this->contig)
		{
			return this->start_pos < ann.start_pos;
		}
		return this->contig < ann.contig;
	}
};

class Overlap
{
public:
	string greatest_overlap;
	int greatest_overlap_length = 0;
	int greatest_overlap_rel_frame = 0;
	vector<string> all_overlaps;
	vector<int> all_overlap_lengths;
	vector<int> all_overlap_rel_frames;
	bool sense_verified_overlap = false;
	bool sense_uncharacterized_overlap = false;
	bool sense_dubious_overlap = false;
	bool sense_pseudogene_overlap = false;
	bool sense_te_overlap = false;
	bool sense_blocked_overlap = false;
	bool trna_overlap = false;
	bool anti_verified_overlap = false;
	bool anti_uncharacterized_overlap = false;
	bool anti_dubious_overlap = false;
	bool anti_pseudogene_overlap = false;
	bool anti_te_overlap = false;
	bool anti_blocked_overlap = false;
	bool ltr_overlap = false;
	bool sno_overlap = false;
	bool ars_overlap = false;
	bool sn_overlap = false;
	bool ncrna_overlap = false;
};

class Ribopos
{
public:
	int chr = -1;
	int pos = -1;
	int dir = -1;
	int nuc = -1;
	int readsum = 0;
	int readlength = 0;
	vector<int> reads;
};

class Experiment
{
public:
	string srr;
	string srp;
	bool preprocessed;
	vector<int> linker;
	map<int, int> frames;
	int paired = 0;
	int end_cutoff = 0;
	int start_cutoff = 0;
	int chx = -1;
	int ypd = -1;
};

struct Pos
{
	int contig = -1;
	int strand = -1;
	int coord = -1;
};

class Ribseq
{
public:
	string id;
	vector<int> seq;
	vector<int> rc_seq;
	vector<int> proc_seq;
	vector<int> phred;
	unsigned long int hash;
	vector<Pos> mappings;
};

class riboseq_scorer
{
public:
	vector<int> hits;
	vector<int> scrambled_hits;
	vector<int> frame_hits;
	double pval = 0;
	double scrambled_pval = 0;
	string best_start_codon = "";
	int best_start_pos = -1;
	double best_start_pval = 0;
	string second_best_start_codon = "";
	int second_best_start_pval = 0;
	double second_best_start_pos = -1;
	riboseq_scorer()
	{
		hits = vector<int>(3);
		scrambled_hits = vector<int>(3);
		frame_hits = vector<int>(3);
	}
};

class MultiBlock
{
public:
	string start_gene;
	string end_gene;
	vector<int> start_pos;
	vector<int> start_gene_start;
	vector<int> start_gene_end;
	vector<int> end_gene_start;
	vector<int> end_gene_end;
	string specific_orf = "X";
	int orf_id = -1;
	vector<int> end_pos;
	vector<int> contig;
	vector<vector<int>> seq;
	vector<vector<int>> align_seq;
	vector<int> in_start_intergene;
	vector<int> in_end_intergene;
	bool single_anchor = false;
	bool align_success = false;
	//vector<Sub> subs;
	//map<int, vector<int>> submap;
	bool is_good = true;
	bool rc = false;
	int blast_id=-1;
	double blast_evalue=1.0;
	MultiBlock(const int SPECIES_COUNT)
	{
		start_pos = vector<int>(SPECIES_COUNT);
		end_pos = vector<int>(SPECIES_COUNT);
		contig = vector<int>(SPECIES_COUNT);
		seq = vector<vector<int>>(SPECIES_COUNT);
		//align_seq = vector<vector<int>>(BRANCH_COUNT);
		align_seq = vector<vector<int>>(SPECIES_COUNT);
		start_gene_start = vector<int>(SPECIES_COUNT);
		start_gene_end = vector<int>(SPECIES_COUNT);
		end_gene_start = vector<int>(SPECIES_COUNT);
		end_gene_end = vector<int>(SPECIES_COUNT);
		in_start_intergene = vector<int>(SPECIES_COUNT);
		in_end_intergene = vector<int>(SPECIES_COUNT);
	}
};

class Vcf
{
public:
	int contig = -1;
	int pos = -1;
	int varcount = 0;
	int max_length = 0;
	int reference_length = 0;
	int major_var = -1;
	//vector<int> reference;
	vector<vector<int>> variants;
	vector<vector<int>> calls;
};

class Hex
{
public:
	int index = 0;
	int trip0 = 0;
	int trip1 = 0;
	int aa0 = -1;
	int aa1 = -1;
	int orf_freq = 0;
	int diaa = 0;
	string aa0_str;
	string aa1_str;
	vector<int> seq;
	vector<int> pre_count;
	vector<int> post_count;
	double score = 0;
	Hex()
	{
		seq = vector<int>(6);
	}
};


struct Blast_info
{
	int id0;
	int id1;
	string str_id0;
	string str_id1;
	double identity = 0.0;
	int align_length = 0;
	int mismatch_count = 0;
	int gap_count = 0;
	int align0_start;
	int align0_end;
	int align1_start;
	int align1_end;
	int strand = -1;
	double evalue = -1;
	double bitscore = -1;
	int length0 = 0;
	int length1 = 0;
};

class BlastMatch
{
public:
	string opp_gene_id;
	string opp_contig_id;
	string contig_id;
	int opp_start_pos=-1;
	int opp_stop_pos = -1;
	int opp_strand = -1;
	vector<int> locus_seq;
	vector<vector<int>> align_seq;
	vector<int> align_pos;
	double eval;
	int length0 = 0;
	int length1 = 0;
	double identity = 0;
	int mismatch_count = 0;;
	int gap_count = 0;
};

class GeneSimilar
{
public:
	string best_match;
	double evalue=100;
	vector<string> all_matches;
	vector<double> all_evalues;
	bool verified_uncharacterized_match = false;
	bool dubious_match = false;
	bool pseudogene_match = false;
	bool te_match = false;
	int identity=0;
};

struct DNDS_reporter
{
public:
	string orf_id = "";
	bool pass = true;
	double syn_exp = 0;
	double nonsyn_exp = 0;
	int syn_obs = 0;
	int nonsyn_obs = 0;
	vector<int> syn_base;
	vector<int> syn_subs;
	vector<int> nonsyn_base;
	vector<int> nonsyn_subs;

	DNDS_reporter()
	{
		syn_base = vector<int>(16);
		syn_subs = vector<int>(16);
		nonsyn_base = vector<int>(16);
		nonsyn_subs = vector<int>(16);
	}
};

struct DNDS_reporter2
{
public:
	string orf_id;
	string context;
	string base;
	string sub;
	int syn = 0;
	int actual = 0;
	string codon;
	string anticodon;
	int codon_pos;;
};

bool is_integer(const string& s)
{
	if (s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false;

	char* p;
	int n = strtol(s.c_str(), &p, 10);

	return (*p == 0);
}

void split(const string& s, char delim, vector<string>& elems)
{
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim)) {
		elems.push_back(item);
	}
}

bool file_exists(string filename)
{
	ifstream file(filename);
	return file.good();
}

void get_filenames(vector<string> &filenames, string filename, bool skip_header)
{
	ifstream file(filename);
	string line;
	if (skip_header)
	{
		getline(file, line);
	}
	while (getline(file, line))
	{
		filenames.push_back(line);
	}
}

void reverse_complement(vector<int>& rc, const vector<int>& seq)
{
	vector<int> comp = { 3,2,1,0,4,5,6 };
	for (int i = seq.size() - 1; i >= 0; i--)
	{
		if (seq[i] > comp.size())
		{
			cout << "\nmismatch:" << seq.size() << " " << seq[i];
			getchar();
		}
		else
		{
			rc.push_back(comp.at(seq[i]));
		}
	}
}

void reverse_complement(vector<int>& seq)
{
	vector<int> rc;
	reverse_complement(rc, seq);
	seq = rc;
}

void read_sgd_annotation(vector<SGD>& sgds, string filename)
{
	const map<string, int> chr_map = { { "chrI",1 },{ "chrII",2 },{ "chrIII",3 },{ "chrIV",4 },{ "chrV",5 },{ "chrVI",6 },{ "chrVII",7 },{ "chrVIII",8 },{ "chrIX",9 },{ "chrX",10 },{ "chrXI",11 },{ "chrXII",12 },{ "chrXIII",13 },{ "chrXIV",14 },{ "chrXV",15 },{ "chrXVI",16 },{ "chrmt",17 } };

	ifstream file(filename);
	string line;
	while (getline(file, line))
	{
		if (line[0] != '#')
		{
			vector<string> columns;
			split(line, '\t', columns);
			sgds.push_back(SGD());
			if (!chr_map.count(columns[0]))
			{
				continue;
			}
			sgds.back().chr = chr_map.at(columns[0]) - 1;
			sgds.back().annotation_type = columns[2];
			if (columns[2] == "pseudogene")
			{
				sgds.back().orf_classification = "pseudogene";
			}
			else if (columns[2] == "transposable_element_gene")
			{
				sgds.back().orf_classification = "transposable_element_gene";
			}
			else if (columns[2] == "LTR_retrotransposon")
			{
				sgds.back().orf_classification = "LTR_retrotransposon";
			}
			else if (columns[2] == "blocked_reading_frame")
			{
				sgds.back().orf_classification = "blocked_reading_frame";
			}
			else if (columns[2] == "tRNA_gene")
			{
				sgds.back().orf_classification = "tRNA_gene";
			}
			else if (columns[2] == "long_terminal_repeat")
			{
				sgds.back().orf_classification = "long_terminal_repeat";
			}
			else if (columns[2] == "snoRNA_gene")
			{
				sgds.back().orf_classification = "snoRNA_gene";
			}
			else if (columns[2] == "snRNA_gene")
			{
				sgds.back().orf_classification = "snRNA_gene";
			}
			else if (columns[2] == "ARS")
			{
				sgds.back().orf_classification = "ARS";
			}
			else if (columns[2] == "ncRNA_gene")
			{
				sgds.back().orf_classification = "ncRNA_gene";
			}

			sgds.back().start = stoi(columns[3]) - 1;
			sgds.back().end = stoi(columns[4]) - 1;
			
			if (columns[6] == "+")
			{
				sgds.back().strand = 0;
			}
			else
			{
				sgds.back().strand = 1;
			}
			vector<string> info;
			split(columns[8], ';', info);
			for (int i = 0; i < info.size(); i++)
			{
				vector<string> components;
				split(info[i], '=', components);
				if (components[0] == "ID")
				{
					sgds.back().ID = components[1];
				}
				if (components[0] == "Name")
				{
					sgds.back().Name = components[1];
				}
				if (components[0] == "Parent")
				{
					sgds.back().Parent = components[1];
				}
				if (components[0] == "orf_classification")
				{
					sgds.back().orf_classification = components[1];
				}
			}
		}
	}
}

void read_fasta(vector<Contig>& contigs, string filename)
{
	ifstream file(filename);
	map<char, int> nucmap = { { 'A',0 },{ 'a',0 },{ 'C',1 },{ 'c',1 },{ 'G',2 },{ 'g',2 },{ 'T',3 },{ 't',3 },{ 'N',5 },{ 'n',5 },{'-',4},{'K',5},{'R',5},{'W',5},{'Y',5},{'k',5},{'r',5},{'w',5},{'y',5},{'S',5},{'s',5} ,{'M',5},{'m',5} ,{'V',5},{'v',5},{'H',5},{'h',5},{'B',5},{'b',5},{'D',5},{'d',5} };
	string line;
	while (getline(file, line))
	{
		if (line.substr(0, 1) == ">")
		{
			contigs.push_back(Contig());
			vector<string> column_data;
			split(line, ' ', column_data);
			contigs.back().id = column_data[0].substr(1);
		}
		else
		{
			for (unsigned int i = 0; i < line.size(); i++)
			{
				if ((int)line[i] != 13)
				{
					if (!nucmap.count(line[i]))
					{
						cout << "\nbad base: " << line << "\n" << line[i] << "\n" << (int)line[i];
						getchar();
					}
					contigs.back().seq.push_back(nucmap.at(line[i]));
				}
			}
		}
	}
}


void read_annotated_genes(vector<Orf>& genes, string filename, bool all_classes)
{
	const map<string, int> roman = { { "I",1 },{ "II",2 },{ "III",3 },{ "IV",4 },{ "V",5 },{ "VI",6 },{ "VII",7 },{ "VIII",8 },{ "IX",9 },{ "X",10 },{ "XI",11 },{ "XII",12 },{ "XIII",13 },{ "XIV",14 },{ "XV",15 },{ "XVI",16 },{ "Mito",17 } };
	ifstream file(filename);
	string line;
	map<char, int> nucmap = { { 'A',0 },{ 'a',0 },{ 'C',1 },{ 'c',1 },{ 'G',2 },{ 'g',2 },{ 'T',3 },{ 't',3 },{ 'N',4 },{ 'n',4 } };
	bool read_pos = false;
	while (getline(file, line))
	{
		if (line.substr(0, 1) == ">")
		{
			read_pos = false;
			string id = "";
			vector<string> column_data;
			split(line, ',', column_data);
			for (unsigned int i = 0; i < column_data.size(); i++)
			{
				vector<string> row_data;
				split(column_data[i], ' ', row_data);
				if (i == 0)
				{
					id = row_data[0].substr(1);
				}
				else if (row_data.size() > 1 && row_data[1] == "Chr")
				{
					string orf_class = "";
					if (line.find(", Verified ORF,") != string::npos)
					{
						orf_class = "verified";
					}
					else if (line.find(", Uncharacterized ORF,") != string::npos)
					{
						orf_class = "uncharacterized";
					}
					else if (line.find(", pseudogene,") != string::npos)
					{
						orf_class = "pseudogene";
					}
					else if (line.find(", telomere,") != string::npos)
					{
						orf_class = "telomere";
					}
					else if (line.find(", blocked_reading_frame,") != string::npos)
					{
						orf_class = "blocked";
					}
					else if ((line.find(", transposable_element_gene,") != string::npos))
					{
						orf_class = "te";
					}
					else if ((line.find(", LTR_retrotransposon,") != string::npos))
					{
						orf_class = "rte";
					}
					else if (line.find(", Dubious ORF,") != string::npos)// Dubious ORF
					{
						orf_class = "dubious";
					}

					if (all_classes || orf_class != "")
					{
						/*if (id == "YNL205C")
						{
							cout << "\nAAA" << id<<" "<<orf_class;
							getchar();
						}*/
						genes.push_back(Orf(1));
						genes.back().orf_class = orf_class;
						genes.back().annotation = id;
						genes.back().annotated = true;
						if (!roman.count(row_data[2]))
						{
							cout << "\nMISS: " << row_data[2];
							getchar();
						}
						genes.back().contig = roman.at(row_data[2]) - 1;
						int pos0 = stoi(row_data[4].substr(0, row_data[4].find('-'))) - 1;
						int pos1 = stoi(row_data[4].substr(row_data[4].find('-') + 1)) - 1;
						if (pos0 > pos1)
						{
							genes.back().strand = 1;
							genes.back().start_pos = pos1;
							genes.back().end_pos = pos0;
						}
						else
						{
							genes.back().strand = 0;
							genes.back().start_pos = pos0;
							genes.back().end_pos = pos1;
						}
						read_pos = true;
					}
				}
				else if (read_pos)
				{
					string pro_num = column_data[i].substr(0, column_data[i].find('-'));
					if (is_integer(pro_num))
					{
						int pos0 = stoi(pro_num) - 1;
						int pos1 = stoi(column_data[i].substr(column_data[i].find('-') + 1)) - 1;
						genes.back().splice = 1;
						if (pos0 > pos1)
						{
							genes.back().end_pos = pos0;
						}
						else
						{
							genes.back().end_pos = pos1;
						}
					}
					else
					{
						break;
					}
				}
			}
		}
		else if (read_pos)
		{
			for (unsigned int i = 0; i < line.size(); i++)
			{
				if (!nucmap.count(line[i]))
				{
					cout << "\nmiss: " << line;
					getchar();
				}
				genes.back().seq.push_back(nucmap.at(line[i]));
			}
		}
	}
}

void organize_sgds(vector<SGD>& sgds_organized, const vector<SGD>& sgds)
{
	map<string, int> gene_to_index;
	map<string, int> mrna_to_index;
	map<string, int> cds_to_index;
	for (int i = 0; i < sgds.size(); i++)
	{
		if (sgds[i].annotation_type == "gene" ||
			sgds[i].annotation_type == "pseudogene" ||
			sgds[i].annotation_type == "transposable_element_gene" ||
			sgds[i].annotation_type == "LTR_retrotransposon" ||
			sgds[i].annotation_type == "tRNA_gene" ||
			sgds[i].annotation_type == "blocked_reading_frame" ||
			sgds[i].annotation_type == "long_terminal_repeat" ||
			sgds[i].annotation_type == "snoRNA_gene" ||
			sgds[i].annotation_type == "snRNA_gene" ||
			sgds[i].annotation_type == "ARS" ||
			sgds[i].annotation_type == "ncRNA_gene"
			)
		{
			sgds_organized.push_back(sgds[i]);
			gene_to_index[sgds[i].Name] = sgds_organized.size() - 1;
		}
	}

	for (int i = 0; i < sgds.size(); i++)
	{
		if (sgds[i].annotation_type == "mRNA")
		{
			if (gene_to_index.count(sgds[i].Parent))
			{
				mrna_to_index[sgds[i].Name] = gene_to_index.at(sgds[i].Parent);
			}
		}
	}

	for (int i = 0; i < sgds.size(); i++)
	{
		if (sgds[i].annotation_type == "CDS")
		{
			if (mrna_to_index.count(sgds[i].Parent))
			{
				int index = mrna_to_index.at(sgds[i].Parent);
				cds_to_index[sgds[i].Name] = index;
				if (sgds_organized[index].cds_start == -1 || sgds_organized[index].cds_start > sgds[i].start)
				{
					sgds_organized[index].cds_start = sgds[i].start;
				}
				if (sgds_organized[index].cds_end == -1 || sgds_organized[index].cds_end < sgds[i].end)
				{
					sgds_organized[index].cds_end = sgds[i].end;
				}
			}
		}
	}

	for (int i = 0; i < sgds.size(); i++)
	{
		if (sgds[i].annotation_type == "intron" || sgds[i].annotation_type == "intein_encoding_region")
		{
			if (mrna_to_index.count(sgds[i].Parent))
			{
				int index = mrna_to_index.at(sgds[i].Parent);
				sgds_organized[index].splice = true;
			}
		}
	}

	for (int i = 0; i < sgds_organized.size(); i++)
	{
		if (sgds_organized[i].cds_start != -1)
		{
			sgds_organized[i].start = sgds_organized[i].cds_start;
		}
		if (sgds_organized[i].cds_end != -1)
		{
			sgds_organized[i].end = sgds_organized[i].cds_end;
		}
	}
}



void rev_comp_contigs(vector<Contig>& contigs)
{
	for (int i = 0; i < contigs.size(); i++)
	{
		vector<int> rc;
		reverse_complement(rc, contigs[i].seq);
		contigs[i].seq = rc;
	}
}

void detect_orfs(vector<Orf>& orfs, const vector<int>& seq, int contig_id, int strand, int species_count)
{
	vector<int> start = { 0,3,2 };
	int offset = 0;
	for (unsigned int j = 0; j < seq.size() - 3; j++)
	{
		if (seq[j] == start[0] && seq[j + 1] == start[1] && seq[j + 2] == start[2])//start codon
		{
			for (int k = j + 3; k < (signed)seq.size() - 3 - offset; k += 3)
			{
				if (seq[k + offset] == 3 && (
					(seq[k + 1 + offset] == 0 && seq[k + 2 + offset] == 0) ||
					(seq[k + 1 + offset] == 2 && seq[k + 2 + offset] == 0) ||
					(seq[k + 1 + offset] == 0 && seq[k + 2 + offset] == 2)
					))
				{
					orfs.push_back(Orf(species_count));
					orfs.back().strand = strand;
					orfs.back().start_pos = j;
					orfs.back().end_pos = k + 2;
					if (orfs.back().strand == 1)
					{
						orfs.back().start_pos = seq.size() - 3 - k;
						orfs.back().end_pos = seq.size() - 1 - (k + 2) + (k + 2 - j);
					}
					orfs.back().contig = contig_id;
					break;
				}
			}
		}
	}
}

void detect_orfs(vector<Orf>& orfs, const vector<Contig>& contigs, int strand, int species_count)
{
	for (unsigned int i = 0; i < contigs.size(); i++)
	{
		detect_orfs(orfs, contigs[i].seq, i, strand, species_count);
	}
}

void filter_orfs_by_length(vector<Orf>& orfs, int thres)
{
	vector<Orf> filtered_orfs;
	for (int i = 0; i < orfs.size(); i++)
	{
		if (orfs[i].end_pos - orfs[i].start_pos >= thres)
		{
			filtered_orfs.push_back(orfs[i]);
		}
	}
	orfs = filtered_orfs;
}

void remove_overlapping_orfs(vector<Orf>& orfs)
{
	vector<int> overlapping(orfs.size());
	for (int i = 0; i < orfs.size() - 1; i++)
	{
		for (int j = i + 1; j < orfs.size(); j++)
		{
			if (orfs[i].end_pos > orfs[j].start_pos&& orfs[i].contig == orfs[j].contig)
			{
				if (abs(orfs[i].end_pos - orfs[i].start_pos) < abs(orfs[j].end_pos - orfs[j].start_pos))
				{
					overlapping[i] = 1;
				}
				else
				{
					overlapping[j] = 1;
				}
			}
			else
			{
				break;
			}
		}
	}
	vector<Orf> filtered_orfs;
	for (int i = 0; i < orfs.size(); i++)
	{
		if (!overlapping[i])
		{
			filtered_orfs.push_back(orfs[i]);
		}
	}
	orfs = filtered_orfs;
}

void get_orf_seqs(vector<Orf>& orfs, const vector<Contig>& contigs)
{
	for (int i = 0; i < orfs.size(); i++)
	{
		orfs[i].seq = vector<int>();
		int contig_id = orfs[i].contig;
		int start = orfs[i].start_pos;
		int end = orfs[i].end_pos;
		for (int j = start; j <= end; j++)
		{
			if ((orfs[i].strand == 0 && j < start + 3) || (orfs[i].strand == 1 && j > end - 3))
			{
				orfs[i].seq.push_back(contigs[contig_id].seq[j]);
			}
			else if (orfs[i].strand == 0)
			{
				orfs[i].seq.push_back(contigs[contig_id].seq[j]);
			}
			else
			{
				orfs[i].seq.push_back(contigs[contig_id].seq[j]);
			}
		}
		if (orfs[i].strand)
		{
			reverse_complement(orfs[i].seq);
		}
	}
}

void remove_orfs_with_n(vector<Orf>& orfs)
{
	vector<Orf> filtered_orfs;
	for (int i = 0; i < orfs.size(); i++)
	{
		bool all_good = true;
		for (int j = 0; j < orfs[i].seq.size(); j++)
		{
			if (orfs[i].seq[j] > 3 || orfs[i].seq[j] < 0)
			{
				all_good = false;
			}
		}
		if (all_good)
		{
			filtered_orfs.push_back(orfs[i]);
		}
	}
	orfs = filtered_orfs;
}

void get_orf_nucseq(vector<Orf>& orfs)
{
	vector<char> nucmap = { 'A','C','G','T','N' };
	for (int i = 0; i < orfs.size(); i++)
	{
		for (int j = 0; j < orfs[i].seq.size(); j++)
		{
			orfs[i].nucseq += nucmap.at(orfs[i].seq[j]);
		}
	}
}

string translate(const vector<int>& seq)
{
	string aaseq;
	map<int, char> aa_map = { { 0,'K' },{ 1,'N' },{ 2,'K' },{ 3,'N' },{ 4,'T' },{ 5,'T' },{ 6,'T' },{ 7,'T' },{ 8,'R' },{ 9,'S' },{ 10,'R' },{ 11,'S' },{ 12,'I' },{ 13,'I' },{ 14,'M' },{ 15,'I' },{ 16,'Q' },{ 17,'H' },{ 18,'Q' },{ 19,'H' },{ 20,'P' },{ 21,'P' },{ 22,'P' },{ 23,'P' },{ 24,'R' },{ 25,'R' },{ 26,'R' },{ 27,'R' },{ 28,'L' },{ 29,'L' },{ 30,'L' },{ 31,'L' },{ 32,'E' },{ 33,'D' },{ 34,'E' },{ 35,'D' },{ 36,'A' },{ 37,'A' },{ 38,'A' },{ 39,'A' },{ 40,'G' },{ 41,'G' },{ 42,'G' },{ 43,'G' },{ 44,'V' },{ 45,'V' },{ 46,'V' },{ 47,'V' },{ 48,'X' },{ 49,'Y' },{ 50,'X' },{ 51,'Y' },{ 52,'S' },{ 53,'S' },{ 54,'S' },{ 55,'S' },{ 56,'X' },{ 57,'C' },{ 58,'W' },{ 59,'C' },{ 60,'L' },{ 61,'F' },{ 62,'L' },{ 63,'F' } };
	for (int i = 0; i < seq.size() - 2; i += 3)
	{
		int hit = seq.at(i) * 4 * 4 + seq.at(i + 1) * 4 + seq.at(i + 2);
		if (!aa_map.count(hit))
		{
			cout << "\ntranslation miss: " << i << " " << seq.size() << " " << seq.at(i) << " " << seq.at(i + 1) << " " << seq.at(i + 2);
			return "";
		}
		aaseq += aa_map.at(hit);
	}
	return aaseq;
}

void translate_genes(vector<Orf>& orfs)
{
	for (int i = 0; i < orfs.size(); i++)
	{
		orfs[i].aaseq = translate(orfs[i].seq);
	}
}

void assign_orf_ids(vector<Orf>& orfs)
{
	for (int i = 0; i < orfs.size(); i++)
	{
		orfs[i].id = i;
		orfs[i].orf_id = to_string(orfs[i].contig) + "_" + to_string(orfs[i].start_pos) + "_" + to_string(orfs[i].end_pos) + "_" + to_string(orfs[i].strand);
	}
}

void get_orfs(vector<Orf>& orfs, const vector<Contig>& contigs, int species_count, int length_filter, bool remove_overlapping)
{
	vector<Contig> rc_contigs_focal = contigs;
	rev_comp_contigs(rc_contigs_focal);

	detect_orfs(orfs, contigs, 0, species_count);
	detect_orfs(orfs, rc_contigs_focal, 1, species_count);
	cout << "\norfs detected: " << orfs.size();
	filter_orfs_by_length(orfs, length_filter);
	cout << "\norfs remaining after size filter: " << orfs.size();

	sort(orfs.begin(), orfs.end());

	//2. eliminate overlapping orfs
	if (remove_overlapping)
	{
		remove_overlapping_orfs(orfs);
		cout << "\norfs remaining after overlap filter: " << orfs.size();
	}

	get_orf_seqs(orfs, contigs);
	remove_orfs_with_n(orfs);
	cout << "\norfs remaining after N filter: " << orfs.size();
	get_orf_nucseq(orfs);
	translate_genes(orfs);
	assign_orf_ids(orfs);
}

void remove_duplicate_orfs_by_contig(vector<Orf>& orfs, vector<Contig>& contigs)
{
	vector<Orf> filtered_orfs;
	vector<map<int, int>> stop_locs(contigs.size());
	for (int i = 0; i < orfs.size(); i++)
	{
		if (orfs[i].strand == 0)
		{
			if (!stop_locs[orfs[i].contig].count(orfs[i].end_pos) || orfs[i].end_pos - orfs[i].start_pos > stop_locs[orfs[i].contig][orfs[i].end_pos])
			{
				stop_locs[orfs[i].contig][orfs[i].end_pos] = orfs[i].end_pos - orfs[i].start_pos;
			}
		}
	}
	for (int i = 0; i < orfs.size(); i++)
	{
		if (orfs[i].strand == 0)
		{
			if (stop_locs[orfs[i].contig][orfs[i].end_pos] == orfs[i].end_pos - orfs[i].start_pos)
			{
				filtered_orfs.push_back(orfs[i]);
			}
		}
	}
	stop_locs = vector<map<int, int>>(contigs.size());
	for (int i = 0; i < orfs.size(); i++)
	{
		if (orfs[i].strand == 1)
		{
			if (!stop_locs[orfs[i].contig].count(orfs[i].start_pos) || orfs[i].end_pos - orfs[i].start_pos > stop_locs[orfs[i].contig][orfs[i].start_pos])
			{
				stop_locs[orfs[i].contig][orfs[i].start_pos] = orfs[i].end_pos - orfs[i].start_pos;
			}
		}
	}
	for (int i = 0; i < orfs.size(); i++)
	{
		if (orfs[i].strand == 1)
		{
			if (stop_locs[orfs[i].contig][orfs[i].start_pos] == orfs[i].end_pos - orfs[i].start_pos)
			{
				filtered_orfs.push_back(orfs[i]);
			}
		}
	}

	orfs = filtered_orfs;
}

int get_max_overlap(int start0, int end0, int start1, int end1)
{
	if (start1 > end0 || start0 > end1)
	{
		return 0;
	}
	else if (start1 >= start0 && end1 <= end0)
	{
		return end1 - start1;
	}
	else if (start1 <= start0 && end1 >= end0)
	{
		return end0 - start0;
	}
	else if (start1 <= start0)
	{
		return end1 - start0;
	}
	else if (start0 <= start1)
	{
		return end0 - start1;
	}
	return 0;
}

void check_orf_overlap(vector<Orf>& orfs, const vector<Orf>& genes)
{
	for (int i = 0; i < orfs.size(); i++)
	{
		for (int j = 0; j < genes.size(); j++)
		{
			if (orfs[i].contig == genes[j].contig &&
				(
				(orfs[i].start_pos <= genes[j].start_pos && orfs[i].end_pos >= genes[j].start_pos) ||
					(orfs[i].start_pos <= genes[j].end_pos && orfs[i].end_pos >= genes[j].end_pos) ||
					(orfs[i].start_pos >= genes[j].start_pos && orfs[i].end_pos <= genes[j].end_pos) ||
					(orfs[i].start_pos <= genes[j].start_pos && orfs[i].end_pos >= genes[j].end_pos)))
			{
				if (orfs[i].start_pos == genes[j].start_pos && orfs[i].end_pos == genes[j].end_pos)
				{
					if (genes[j].annotated)
					{
						orfs[i].is_gene = genes[j].annotation;
						orfs[i].orf_class = genes[j].orf_class;
					}
					else
					{
						orfs[i].is_gene = to_string(genes[j].id);
					}
				}
				else if (genes[j].splice)
				{
					orfs[i].splice_gene = genes[j].annotation;
				}
				else if ((orfs[i].strand == 0 && orfs[i].end_pos == genes[j].end_pos && orfs[i].start_pos < genes[j].start_pos) ||
					(orfs[i].strand == 1 && orfs[i].start_pos == genes[j].start_pos && orfs[i].end_pos > genes[j].end_pos))
				{
					orfs[i].up_extend_gene = genes[j].annotation;

				}
				else
				{
					if (genes[j].annotated)
					{
						orfs[i].overlaps_gene = genes[j].annotation;
						orfs[i].overlap_gene_relative_frame = 1 + (abs(orfs[i].start_pos - genes[j].start_pos) % 3);
						if (orfs[i].strand != genes[j].strand)
						{
							orfs[i].overlap_gene_relative_frame *= -1;
						}
						orfs[i].overlap_gene_length = get_max_overlap(orfs[i].start_pos, orfs[i].end_pos, genes[j].start_pos, genes[j].end_pos);
					}
					else
					{
						orfs[i].overlaps_gene = to_string(genes[j].id);
					}
				}
			}
		}
	}
}

void check_orf_overlap(vector<Orf>& orfs, const vector<SGD>& sgds, vector<Overlap>& overlaps)
{
	overlaps = vector<Overlap>(orfs.size());
	int start_index = 0;
	for (int i = 0; i < orfs.size(); i++)
	{
		while (sgds[start_index].chr < orfs[i].contig || orfs[i].start_pos - sgds[start_index].end>1000000)
		{
			start_index++;
		}
		for (int j = start_index; j < sgds.size(); j++)
		{
			if (sgds[j].chr > orfs[i].contig || sgds[j].start > orfs[i].end_pos)
			{
				break;
			}
			if (orfs[i].contig == sgds[j].chr &&
				(
				(orfs[i].start_pos <= sgds[j].start && orfs[i].end_pos >= sgds[j].start) ||
					(orfs[i].start_pos <= sgds[j].end && orfs[i].end_pos >= sgds[j].end) ||
					(orfs[i].start_pos >= sgds[j].start && orfs[i].end_pos <= sgds[j].end) ||
					(orfs[i].start_pos <= sgds[j].start && orfs[i].end_pos >= sgds[j].end)))
			{
				if (orfs[i].start_pos == sgds[j].start && orfs[i].end_pos == sgds[j].end)
				{
					orfs[i].is_gene = sgds[j].Name;
					orfs[i].orf_class = sgds[j].orf_classification;
				}
				else if (sgds[j].splice)
				{
					orfs[i].splice_gene = sgds[j].Name;
				}
				else if ((orfs[i].strand == 0 && orfs[i].end_pos == sgds[j].end && orfs[i].start_pos < sgds[j].start) ||
					(orfs[i].strand == 1 && orfs[i].start_pos == sgds[j].start && orfs[i].end_pos > sgds[j].end))
				{
					orfs[i].up_extend_gene = sgds[j].Name;
					if (orfs[i].orf_class == "None")
					{
						orfs[i].orf_class = sgds[j].orf_classification;
					}
				}
				else //overlap
				{
					orfs[i].overlaps.push_back(OverlapHit());
					orfs[i].overlaps.back().annotation = sgds[j].Name;
					orfs[i].overlaps.back().classification = sgds[j].orf_classification;
					orfs[i].overlaps.back().start_pos = sgds[j].start;
					orfs[i].overlaps.back().end_pos = sgds[j].end;
					if ((orfs[i].strand == 0 && orfs[i].start_pos > sgds[j].start)||
						(orfs[i].strand == 1 && orfs[i].end_pos < sgds[j].end ))
					{
						orfs[i].start_codon_overlap = true;
					}
					if ((orfs[i].strand == 0 && orfs[i].end_pos < sgds[j].end) ||
						(orfs[i].strand == 1 && orfs[i].start_pos > sgds[j].start))
					{
						orfs[i].stop_codon_overlap = true;
					}

					//orfs[i].overlaps_gene = sgds[j].Name;

					int rel_frame =  1+((3000000 + orfs[i].start_pos - sgds[j].start) % 3);
					if (rel_frame == 1)
					{
						rel_frame = 3;
					}
					else if (rel_frame == 3)
					{
						rel_frame = 1;
					}
					if (orfs[i].strand != sgds[j].strand)
					{
						rel_frame *= -1;
					}

					/*int rel_frame = 1 + (abs(orfs[i].start_pos - sgds[j].start) % 3);
					if (orfs[i].strand != sgds[j].strand)
					{
						rel_frame *= -1;
					}*/
					//orfs[i].overlap_gene_relative_frame = rel_frame;
					int overlap_length = get_max_overlap(orfs[i].start_pos, orfs[i].end_pos, sgds[j].start, sgds[j].end);
					//cout << "\nol:" << i << " " << overlap_length << " " << orfs[i].start_pos << " " << orfs[i].end_pos << " " << sgds[j].start << " " << sgds[j].end;
					//getchar();
					//orfs[i].overlap_gene_length = overlap_length;
					if (rel_frame<0 && overlap_length > orfs[i].overlap_gene_length)
					{
						orfs[i].overlap_gene_relative_frame = rel_frame;
						orfs[i].overlaps_gene = sgds[j].Name;
						orfs[i].overlap_gene_length = overlap_length;
					}

					if (overlap_length > overlaps[i].greatest_overlap_length)
					{
						overlaps[i].greatest_overlap_length = overlap_length;
						overlaps[i].greatest_overlap = sgds[j].Name;
						overlaps[i].greatest_overlap_rel_frame = rel_frame;
					}
					overlaps[i].all_overlaps.push_back(sgds[j].Name);
					overlaps[i].all_overlap_lengths.push_back(overlap_length);
					overlaps[i].all_overlap_rel_frames.push_back(rel_frame);
					if (sgds[j].orf_classification == "Verified" && rel_frame > 0)
					{
						overlaps[i].sense_verified_overlap = 1;
					}
					if (sgds[j].orf_classification == "Verified" && rel_frame < 0)
					{
						overlaps[i].anti_verified_overlap = 1;
					}
					if (sgds[j].orf_classification == "Uncharacterized" && rel_frame > 0)
					{
						overlaps[i].sense_uncharacterized_overlap = 1;
					}
					if (sgds[j].orf_classification == "Uncharacterized" && rel_frame < 0)
					{
						overlaps[i].anti_uncharacterized_overlap = 1;
					}
					if (sgds[j].orf_classification == "Dubious" && rel_frame > 0)
					{
						overlaps[i].sense_dubious_overlap = 1;
					}
					if (sgds[j].orf_classification == "Dubious" && rel_frame < 0)
					{
						overlaps[i].anti_dubious_overlap = 1;
					}
					if (sgds[j].orf_classification == "pseudogene" && rel_frame > 0)
					{
						overlaps[i].sense_pseudogene_overlap = 1;
					}
					if (sgds[j].orf_classification == "pseudogene" && rel_frame < 0)
					{
						overlaps[i].anti_pseudogene_overlap = 1;
					}
					if (sgds[j].orf_classification == "transposable_element_gene" && rel_frame > 0)
					{
						overlaps[i].sense_te_overlap = 1;
					}
					if (sgds[j].orf_classification == "transposable_element_gene" && rel_frame < 0)
					{
						overlaps[i].anti_te_overlap = 1;
					}
					if (sgds[j].orf_classification == "LTR_retrotransposon" && rel_frame > 0)
					{
						overlaps[i].sense_te_overlap = 1;
					}
					if (sgds[j].orf_classification == "LTR_retrotransposon" && rel_frame < 0)
					{
						overlaps[i].anti_te_overlap = 1;
					}
					if (sgds[j].orf_classification == "blocked_reading_frame" && rel_frame > 0)
					{
						overlaps[i].sense_blocked_overlap = 1;
					}
					if (sgds[j].orf_classification == "blocked_reading_frame" && rel_frame < 0)
					{
						overlaps[i].anti_blocked_overlap = 1;
					}
					if (sgds[j].orf_classification == "tRNA_gene")
					{
						overlaps[i].trna_overlap = 1;
					}
					if (sgds[j].orf_classification == "long_terminal_repeat")
					{
						overlaps[i].ltr_overlap = 1;
					}
					if (sgds[j].orf_classification == "snoRNA_gene")
					{
						overlaps[i].sno_overlap = 1;
					}
					if (sgds[j].orf_classification == "snRNA_gene")
					{
						overlaps[i].sn_overlap = 1;
					}
					if (sgds[j].orf_classification == "ARS")
					{
						overlaps[i].ars_overlap = 1;
					}
					if (sgds[j].orf_classification == "ncRNA_gene")
					{
						overlaps[i].ncrna_overlap = 1;
					}
				}
			}
		}
	}
}

void read_tifseqs(vector<Tifseq>& tseqs, string filename)
{
	ifstream file(filename);
	string line;
	while (getline(file, line))
	{
		tseqs.push_back(Tifseq());
		vector<string> columns;
		split(line, '\t', columns);
		if (is_integer(columns[0]))
		{
			tseqs.back().chr = stoi(columns[0]) - 1;
			tseqs.back().start = stoi(columns[2]);
			tseqs.back().end = stoi(columns[3]);
			tseqs.back().strand = 0;
			if (tseqs.back().start > tseqs.back().end)
			{
				tseqs.back().start = stoi(columns[3]);
				tseqs.back().end = stoi(columns[2]);
				tseqs.back().strand = 1;
			}
			tseqs.back().hits0 = stoi(columns[4]);
			tseqs.back().hits1 = stoi(columns[5]);
		}
	}
}

void get_tifseq_context(vector<SharedTranscript>& shared, vector<Orf>& orfs, vector<SGD>& sgds, vector<Tifseq>& tseqs)
{
	set<string> relations;
	for (int i = 0; i < tseqs.size(); i++)
	{
		for (int j = 0; j < tseqs[i].covered_genes.size(); j++)
		{
			for (int k = 0; k < tseqs[i].covered_orfs.size(); k++)
			{
				int orf = tseqs[i].covered_orfs[k];
				int gene = tseqs[i].covered_genes[j];
				orfs[orf].coupled_genes.push_back(sgds[gene].Name);
				string orientation;
				if ((orfs[orf].start_pos < sgds[gene].start && sgds[gene].strand == 0) || (orfs[orf].end_pos > sgds[gene].end&& sgds[gene].strand == 1))
				{
					orientation = "upstream";
				}
				else if ((orfs[orf].start_pos < sgds[gene].start && sgds[gene].strand == 1) || (orfs[orf].end_pos > sgds[gene].end&& sgds[gene].strand == 0))
				{
					orientation = "downstream";
				}
				else
				{
					orientation = "within";
				}
				orfs[orf].transcript_relation = orientation;
				string relation = orientation + orfs[orf].orf_id + "_" + sgds[gene].Name;
				if (!relations.count(relation))
				{
					relations.insert(relation);
					shared.push_back(SharedTranscript());
					shared.back().orf_id = orfs[orf].orf_id;
					shared.back().relation = orientation;
					shared.back().shared_sgd = sgds[gene].Name;
					shared.back().shared_sgd_class = sgds[gene].annotation_type;
				}

			}
		}
	}
}

void get_tifseq_coverage(vector<Orf>& loci, vector<Tifseq>& tseqs, bool genes)
{
	int start_index = 0;
	for (int i = 0; i < loci.size(); i++)
	{
		while (tseqs[start_index].chr < loci[i].contig || loci[i].start_pos - tseqs[start_index].end>1000000)
		{
			if (start_index >= tseqs.size())
			{
				break;
			}
			start_index++;
		}
		for (int j = start_index; j < tseqs.size(); j++)
		{
			if (tseqs[j].chr > loci[i].contig || tseqs[j].start > loci[i].end_pos)
			{
				break;
			}
			if (loci[i].contig == tseqs[j].chr && loci[i].strand == tseqs[j].strand && loci[i].start_pos >= tseqs[j].start && loci[i].end_pos <= tseqs[j].end)
			{
				loci[i].tifseqs_hits0 += tseqs[j].hits0;
				loci[i].tifseqs_hits1 += tseqs[j].hits1;
				if (genes)
				{
					tseqs[j].covered_genes.push_back(i);
				}
				else
				{
					tseqs[j].covered_orfs.push_back(i);
				}
			}
		}
	}
}

void get_tifseq_coverage(vector<SGD>& loci, vector<Tifseq>& tseqs, bool genes)
{
	int start_index = 0;
	for (int i = 0; i < loci.size(); i++)
	{
		//cout << "\nloci: " << i;
		while (tseqs[start_index].chr < loci[i].chr || loci[i].start - tseqs[start_index].end>1000000)
		{
			if (start_index >= tseqs.size())
			{
				break;
			}
			start_index++;
		}
		for (int j = start_index; j < tseqs.size(); j++)
		{
			if (tseqs[j].chr > loci[i].chr || tseqs[j].start > loci[i].end)
			{
				break;
			}
			if (loci[i].chr == tseqs[j].chr && loci[i].strand == tseqs[j].strand && loci[i].start >= tseqs[j].start && loci[i].end <= tseqs[j].end)
			{
				loci[i].tifseqs_hits0 += tseqs[j].hits0;
				loci[i].tifseqs_hits1 += tseqs[j].hits1;
				if (genes)
				{
					tseqs[j].covered_genes.push_back(i);
				}
				else
				{
					tseqs[j].covered_orfs.push_back(i);
				}
			}
		}
	}
}

void sort_orfs_by_strand(vector<Orf>& orfs, vector<Overlap>& overlaps)
{
	vector<Orf> rearranged_orfs;
	vector<Overlap> rearranged_overlaps;
	for (int i = 0; i < orfs.size(); i++)
	{
		if (orfs[i].strand == 0)
		{
			rearranged_orfs.push_back(orfs[i]);
			rearranged_overlaps.push_back(overlaps[i]);
		}
	}
	for (int i = 0; i < orfs.size(); i++)
	{
		if (orfs[i].strand == 1)
		{
			rearranged_orfs.push_back(orfs[i]);
			rearranged_overlaps.push_back(overlaps[i]);
		}
	}
	orfs = rearranged_orfs;
	overlaps = rearranged_overlaps;
}

void print_orfs_aa_fasta(string filename, const vector<Orf>& orfs)
{
	ofstream file(filename);
	for (int i = 0; i < orfs.size(); i++)
	{
		file << ">" << "A" << orfs[i].id << "\n";
		file << orfs[i].aaseq.substr(0, orfs[i].aaseq.size()) << "\n";
	}
}

void print_orfs_aa_fasta(string filename, const vector<Orf>& orfs, bool use_orf_id)
{
	ofstream file(filename);
	for (int i = 0; i < orfs.size(); i++)
	{
		if (use_orf_id)
		{
			file << ">" << orfs[i].orf_id << "\n";
		}
		else
		{
			file << ">seq" << i << "\n";
		}
		file << orfs[i].aaseq.substr(0, orfs[i].aaseq.size() - 1) << "\n";
	}
}

void print_orfs_aa_fasta_scrambled_and_real(string filename, const vector<Orf>& orfs, bool use_orf_id)
{
	ofstream file(filename);
	for (int i = 0; i < orfs.size(); i++)
	{
		if (use_orf_id)
		{
			file << ">" << orfs[i].orf_id << "\n";
		}
		else
		{
			file << ">seq" << i << "\n";
		}
		file << orfs[i].aaseq.substr(0, orfs[i].aaseq.size() - 1) << "\n";

		if (use_orf_id)
		{
			file << ">" << orfs[i].orf_id;
		}
		else
		{
			file << ">seq" << i;
		}
		file << "_X\n";
		string aaseq_scrambled = orfs[i].aaseq.substr(0, orfs[i].aaseq.size() - 1);
		random_shuffle(aaseq_scrambled.begin() + 1, aaseq_scrambled.end());
		file << aaseq_scrambled << "\n";
	}
}

void print_orfs_nuc_fasta(string filename, const vector<Orf>& orfs)
{
	ofstream file(filename);
	for (int i = 0; i < orfs.size(); i++)
	{
		file << ">" << orfs[i].orf_id << "\n";
		file << orfs[i].nucseq << "\n";
	}
}
void print_orfs(const vector<Orf>& orfs, string suffix)
{
	ofstream file("orfs" + suffix);
	file << "id orf_id contig start end strand overlaps_gene gene_overlap_length gene_rel_frame is_gene splice_gene up_extend_gene tseqs0 tseqs1 orf_class";
	for (int i = 0; i < orfs.size(); i++)
	{
		file << "\n" << i << " " << orfs[i].orf_id<<" "<< orfs[i].contig << " " << orfs[i].start_pos << " " << orfs[i].end_pos << " " << orfs[i].strand;
		file << " " << orfs[i].overlaps_gene << " " << orfs[i].overlap_gene_length << " " << orfs[i].overlap_gene_relative_frame;
		file << " " << orfs[i].is_gene << " " << orfs[i].splice_gene << " " << orfs[i].up_extend_gene;
		file << " " << orfs[i].tifseqs_hits0 << " " << orfs[i].tifseqs_hits1<<" "<<orfs[i].orf_class; 
	}
}


void print_matched_orfs(const vector<Orf>& orfs, string suffix)
{
	ofstream file("orfs" + suffix);
	map<int, char> rev_nuc_map = { { 0,'A' },{ 1,'C' },{ 2,'G' },{ 3,'T' },{ 4,'-' },{ 5,'N' },{ 6,'N' } };
	map<int, char> match_frame_map = { {4,'-'},{0,'N'},{1,'Y'} };
	file << "id contig start end strand block_id within_block_id orf_id";
	file << " syn_vars nonsyn_vars syn_exp nonsyn_exp nuc_diverse";
	for (int i = 0; i < orfs[0].good_align.size(); i++)
	{
		file << " good" << i;
	}
	for (int i = 0; i < orfs[0].overlap.size(); i++)
	{
		file << " overlap" << i;
	}
	for (int i = 0; i < orfs[0].blast_align_length.size(); i++)
	{
		file << " blast_align_length" << i;
	}
	for (int i = 0; i < orfs[0].blast_evalue.size(); i++)
	{
		file << " evalue" << i;
	}
	for (int i = 0; i < orfs[0].syn_diffs.size(); i++)
	{
		file << " syn_diffs" << i;
	}
	for (int i = 0; i < orfs[0].nonsyn_diffs.size(); i++)
	{
		file << " nonsyn_diffs" << i;
	}
	for (int i = 0; i < orfs[0].syn_diffs_exp.size(); i++)
	{
		file << " syn_diffs_exp" << i;
	}
	for (int i = 0; i < orfs[0].nonsyn_diffs_exp.size(); i++)
	{
		file << " nonsyn_diffs_exp" << i;
	}
	for (int i = 0; i < orfs[0].homo_length.size(); i++)
	{
		file << " homo_length" << i;
	}
	for (int i = 0; i < orfs[0].focal_align_seq.size(); i++)
	{
		file << " focal_align" << i;
	}
	for (int i = 0; i < orfs[0].align_seq.size(); i++)
	{
		file << " align" << i;
	}
	for (int i = 0; i < orfs[0].match_frame.size(); i++)
	{
		file << " match_frame" << i;
	}
	for (int i = 0; i < orfs[0].align_identity.size(); i++)
	{
		file << " align_identity" << i;
	}
	//file << " out0 out1";
	for (int i = 0; i < orfs[0].align_evalues.size(); i++)
	{
		file << " align_evalues" << i;
	}
	file << " focal_aaseq";
	for (int i = 0; i < orfs[0].alt_aaseq.size(); i++)
	{
		file << " aaseq" << i;
	}
	for (int i = 0; i < orfs[0].alt_aaseq.size(); i++)
	{
		file << " focal_aaseq_align" << i;
	}
	for (int i = 0; i < orfs[0].alt_aaseq.size(); i++)
	{
		file << " alt_aaseq_align" << i;
	}
	file << " focal_seq";
	for (int i = 0; i < orfs[0].alt_seqs.size(); i++)
	{
		file << " seq" << i;
	}
	file << " repeats orf_class";
	for (int i = 0; i < orfs[0].focal_flank_align.size(); i++)
	{
		file << " focal_flank" << i;
	}
	for (int i = 0; i < orfs[0].flank_align.size(); i++)
	{
		file << " flank" << i;
	}
	for (int i = 0; i < orfs[0].use_blast_align.size(); i++)
	{
		file << " use_blast" << i;
	}
	for (int i = 0; i < orfs[0].homology_validated.size(); i++)
	{
		file << " homology_validated" << i;
	}
	for (int i = 0; i < orfs.size(); i++)
	{
		file << "\n" << i << " " << orfs[i].contig << " " << orfs[i].start_pos << " " << orfs[i].end_pos << " " << orfs[i].strand;
		file << " " << orfs[i].block_id << " " << orfs[i].within_block_id << " "<<orfs[i].orf_id;
		file << " " << orfs[i].syn_vars << " " << orfs[i].nonsyn_vars << " " << orfs[i].syn_exp << " " << orfs[i].nonsyn_exp << " " << orfs[i].nuc_diverse;
		for (int j = 0; j < orfs[i].good_align.size(); j++)
		{
			file << " " << orfs[i].good_align[j];
		}
		for (int j = 0; j < orfs[i].overlap.size(); j++)
		{
			file << " " << orfs[i].overlap[j];
		}
		for (int j = 0; j < orfs[i].blast_align_length.size(); j++)
		{
			file << " " << orfs[i].blast_align_length[j];
		}
		for (int j = 0; j < orfs[i].blast_evalue.size(); j++)
		{
			file << " " << orfs[i].blast_evalue[j];
		}
		for (int j = 0; j < orfs[i].syn_diffs.size(); j++)
		{
			file << " " << orfs[i].syn_diffs[j];
		}
		for (int j = 0; j < orfs[i].nonsyn_diffs.size(); j++)
		{
			file << " " << orfs[i].nonsyn_diffs[j];
		}
		for (int j = 0; j < orfs[i].syn_diffs_exp.size(); j++)
		{
			file << " " << orfs[i].syn_diffs_exp[j];
		}
		for (int j = 0; j < orfs[i].nonsyn_diffs_exp.size(); j++)
		{
			file << " " << orfs[i].nonsyn_diffs_exp[j];
		}
		for (int j = 0; j < orfs[i].homo_length.size(); j++)
		{
			file << " " << orfs[i].homo_length[j];
		}
		for (int j = 0; j < orfs[i].focal_align_seq.size(); j++)
		{
			file << " ";
			for (int k = 0; k < orfs[i].focal_align_seq[j].size(); k++)
			{
				file << rev_nuc_map.at(orfs[i].focal_align_seq[j][k]);
			}
		}

		for (int j = 0; j < orfs[i].align_seq.size(); j++)
		{
			file << " ";
			for (int k = 0; k < orfs[i].align_seq[j].size(); k++)
			{
				file << rev_nuc_map.at(orfs[i].align_seq[j][k]);
			}
		}
		for (int j = 0; j < orfs[i].match_frame.size(); j++)
		{
			file << " ";
			for (int k = 0; k < orfs[i].match_frame[j].size(); k++)
			{
				file << match_frame_map.at(orfs[i].match_frame[j][k]);
			}
		}
		for (int j = 0; j < orfs[i].align_identity.size(); j++)
		{
			file << " " << orfs[i].align_identity[j];
		}
		//file << " " << orfs[i].outgroup[0] << " " << orfs[i].outgroup[1];
		for (int j = 0; j < orfs[i].align_evalues.size(); j++)
		{
			file << " " << orfs[i].align_evalues[j];
		}
		file << " " << orfs[i].aaseq;
		for (int j = 0; j < orfs[i].alt_aaseq.size(); j++)
		{
			file << " " << orfs[i].alt_aaseq[j];
		}
		for (int j = 0; j < orfs[i].alt_aaseq.size(); j++)
		{
			file << " " << orfs[i].focal_align_aaseq[j];
		}
		for (int j = 0; j < orfs[i].alt_aaseq.size(); j++)
		{
			file << " " << orfs[i].alt_align_aaseq[j];
		}
		file << " ";
		for (int j = 0; j < orfs[i].seq.size(); j++)
		{
			file << rev_nuc_map.at(orfs[i].seq[j]);
		}
		for (int j = 0; j < orfs[i].alt_seqs.size(); j++)
		{
			file << " ";
			for (int k = 0; k < orfs[i].alt_seqs[j].size(); k++)
			{
				file << rev_nuc_map.at(orfs[i].alt_seqs[j][k]);
			}
		}
		file << " " << orfs[i].repeat_count << " " << orfs[i].orf_class;
		for (int j = 0; j < orfs[i].focal_flank_align.size(); j++)
		{
			file << " ";
			for (int k = 0; k < orfs[i].focal_flank_align[j].size(); k++)
			{
				if (!rev_nuc_map.count(orfs[i].focal_flank_align[j][k]))
				{
					cout << "\nhuha: " << i << " " << j << " " << k << " " << orfs[i].focal_flank_align[j][k] << "\ntrial:" << orfs[236097].focal_flank_align[1][8646];
					getchar();
				}
				else
				{
					file << rev_nuc_map.at(orfs[i].focal_flank_align[j][k]);
				}
			}
		}
		for (int j = 0; j < orfs[i].flank_align.size(); j++)
		{
			file << " ";
			for (int k = 0; k < orfs[i].flank_align[j].size(); k++)
			{
				if (!rev_nuc_map.count(orfs[i].flank_align[j][k]))
				{
				}
				else
				{
					file << rev_nuc_map.at(orfs[i].flank_align[j][k]);
				}
			}
		}
		//file << " " << orfs[i].start_codon_overlap << " " << orfs[i].stop_codon_overlap;
		for (int j = 0; j < orfs[i].use_blast_align.size(); j++)
		{
			file << " " << orfs[i].use_blast_align[j];
		}
		for (int j = 0; j < orfs[i].homology_validated.size(); j++)
		{
			file << " "<<orfs[i].homology_validated[j];
		}		
	}
}

void read_matched_orfs(vector<Orf>& orfs, string filename, int species_num)
{
	map<int, char> nuc_map = { { 'A',0 },{ 'C',1 },{ 'G',2 },{ 'T',3 },{ '-',4 },{ 'N',5 } };
	map<char, int> match_frame_map = { {'-',4},{'N',0},{'Y',1} };
	ifstream file(filename);
	string line;
	getline(file, line);
	while (getline(file, line))
	{
		int i = 1;
		vector<string> cols;
		split(line, ' ', cols);
		orfs.push_back(Orf(species_num));
		orfs.back().contig = stoi(cols[i++]);
		orfs.back().start_pos = stoi(cols[i++]);
		orfs.back().end_pos = stoi(cols[i++]);
		orfs.back().strand = stoi(cols[i++]);
		orfs.back().block_id = stoi(cols[i++]);
		orfs.back().within_block_id = stoi(cols[i++]);
		orfs.back().orf_id = cols[i++];

		orfs.back().syn_vars = stoi(cols[i++]);
		orfs.back().nonsyn_vars = stoi(cols[i++]);
		orfs.back().syn_exp = stod(cols[i++]);
		orfs.back().nonsyn_exp = stod(cols[i++]);
		orfs.back().nuc_diverse = stod(cols[i++]);

		for (int j = 0; j < species_num; j++)
		{
			orfs.back().good_align[j] = stoi(cols[i++]);
		}
		for (int j = 0; j < species_num; j++)
		{
			orfs.back().overlap[j] = stoi(cols[i++]);
		}
		for (int j = 0; j < species_num; j++)
		{
			orfs.back().blast_align_length[j] = stoi(cols[i++]);
		}
		for (int j = 0; j < species_num; j++)
		{
			orfs.back().blast_evalue[j] = stod(cols[i++]);
		}

		for (int j = 0; j < species_num; j++)
		{
			orfs.back().syn_diffs[j] = stoi(cols[i++]);
		}
		for (int j = 0; j < species_num; j++)
		{
			orfs.back().nonsyn_diffs[j] = stoi(cols[i++]);
		}
		for (int j = 0; j < species_num; j++)
		{
			orfs.back().syn_diffs_exp[j] = stod(cols[i++]);
		}
		for (int j = 0; j < species_num; j++)
		{
			orfs.back().nonsyn_diffs_exp[j] = stod(cols[i++]);
		}

		for (int j = 0; j < species_num; j++)
		{
			orfs.back().homo_length[j] = stoi(cols[i++]);
		}

		for (int j = 0; j < species_num; j++)
		{
			for (int k = 0; k < cols[i].size(); k++)
			{
				orfs.back().focal_align_seq[j].push_back(nuc_map.at(cols[i][k]));
			}
			i++;
		}

		for (int j = 0; j < species_num; j++)
		{
			for (int k = 0; k < cols[i].size(); k++)
			{
				orfs.back().align_seq[j].push_back(nuc_map.at(cols[i][k]));
			}
			i++;
		}

		for (int j = 0; j < species_num && i < cols.size(); j++)
		{
			if (cols[i] != "NA")
			{
				for (int k = 0; k < cols[i].size(); k++)
				{
					orfs.back().match_frame[j].push_back(match_frame_map.at(cols[i][k]));
				}
			}
			i++;
		}
		for (int j = 0; j < species_num && i < cols.size(); j++)
		{
			orfs.back().align_identity[j] = stod(cols[i]);
			i++;
		}
		for (int j = 0; j < species_num && i < cols.size(); j++)
		{
			orfs.back().align_evalues[j] = stod(cols[i]);
			i++;
		}
		orfs.back().aaseq = cols[i];
		i++;
		for (int j = 0; j < species_num && i < cols.size(); j++)
		{
			orfs.back().alt_aaseq[j] = cols[i];
			i++;
		}
		for (int j = 0; j < species_num && i < cols.size(); j++)
		{
			orfs.back().focal_align_aaseq[j] = cols[i];
			i++;
		}
		for (int j = 0; j < species_num && i < cols.size(); j++)
		{
			orfs.back().alt_align_aaseq[j] = cols[i];
			i++;
		}
		for (int k = 0; k < cols[i].size(); k++)
		{
			orfs.back().seq.push_back(nuc_map.at(cols[i][k]));
		}
		i++;
		for (int j = 0; j < species_num && i < cols.size(); j++)
		{
			for (int k = 0; k < cols[i].size(); k++)
			{
				orfs.back().alt_seqs[j].push_back(nuc_map.at(cols[i][k]));
			}
			i++;
		}
		orfs.back().repeat_count = stoi(cols[i++]);
		orfs.back().orf_class = cols[i++];
		for (int j = 0; j < species_num && i < cols.size(); j++)
		{
			for (int k = 0; k < cols[i].size(); k++)
			{
				orfs.back().focal_flank_align[j].push_back(nuc_map.at(cols[i][k]));
			}
			i++;
		}
		for (int j = 0; j < species_num && i < cols.size(); j++)
		{
			for (int k = 0; k < cols[i].size(); k++)
			{
				orfs.back().flank_align[j].push_back(nuc_map.at(cols[i][k]));
			}
			i++;
		}
		if (i >= cols.size())
		{
			continue;
		}
		for (int j = 0; j < species_num; j++)
		{
			orfs.back().use_blast_align[j] = stoi(cols[i++]);
		}
	}
}


void read_orfs_long(vector<Orf>& orfs, string filename, int species_num)
{
	map<int, char> nuc_map = { { 'A',0 },{ 'C',1 },{ 'G',2 },{ 'T',3 },{ '-',4 },{ 'N',5 } };
	map<char, int> match_frame_map = { {'-',4},{'N',0},{'Y',1} };
	ifstream file(filename);
	string line;
	getline(file, line);
	while (getline(file, line))
	{
		int i = 1;
		vector<string> cols;
		split(line, ' ', cols);
		orfs.push_back(Orf(species_num));
		orfs.back().contig = stoi(cols[i++]);
		orfs.back().start_pos = stoi(cols[i++]);
		orfs.back().end_pos = stoi(cols[i++]);
		orfs.back().strand = stoi(cols[i++]);
		orfs.back().block_id = stoi(cols[i++]);
		orfs.back().within_block_id = stoi(cols[i++]);
		orfs.back().orf_id = cols[i++];
		orfs.back().syn_vars = stoi(cols[i++]);
		orfs.back().nonsyn_vars = stoi(cols[i++]);
		orfs.back().syn_exp = stod(cols[i++]);
		orfs.back().nonsyn_exp = stod(cols[i++]);
		orfs.back().nuc_diverse = stod(cols[i++]);

		for (int j = 0; j < species_num; j++)
		{
			orfs.back().good_align[j] = stoi(cols[i++]);
		}
		for (int j = 0; j < species_num; j++)
		{
			orfs.back().overlap[j] = stoi(cols[i++]);
		}
		for (int j = 0; j < species_num; j++)
		{
			orfs.back().blast_align_length[j] = stoi(cols[i++]);
		}
		for (int j = 0; j < species_num; j++)
		{
			orfs.back().blast_evalue[j] = stod(cols[i++]);
		}

		for (int j = 0; j < species_num; j++)
		{
			orfs.back().syn_diffs[j] = stoi(cols[i++]);
		}
		for (int j = 0; j < species_num; j++)
		{
			orfs.back().nonsyn_diffs[j] = stoi(cols[i++]);
		}
		for (int j = 0; j < species_num; j++)
		{
			orfs.back().syn_diffs_exp[j] = stod(cols[i++]);
		}
		for (int j = 0; j < species_num; j++)
		{
			orfs.back().nonsyn_diffs_exp[j] = stod(cols[i++]);
		}

		for (int j = 0; j < species_num; j++)
		{
			orfs.back().homo_length[j] = stoi(cols[i++]);
		}

		for (int j = 0; j < species_num; j++)
		{
			for (int k = 0; k < cols[i].size(); k++)
			{
				orfs.back().focal_align_seq[j].push_back(nuc_map.at(cols[i][k]));
			}
			i++;
		}
		for (int j = 0; j < species_num; j++)
		{
			for (int k = 0; k < cols[i].size(); k++)
			{
				orfs.back().align_seq[j].push_back(nuc_map.at(cols[i][k]));
			}
			i++;
		}

		for (int j = 0; j < species_num && i<cols.size(); j++)
		{

			if (cols[i] != "NA")
			{
				for (int k = 0; k < cols[i].size(); k++)
				{
					orfs.back().match_frame[j].push_back(match_frame_map.at(cols[i][k]));
				}
			}
			i++;
		}
		for (int j = 0; j < species_num && i < cols.size(); j++)
		{
			orfs.back().align_identity[j] = stod(cols[i]);
			i++;
		}
		for (int j = 0; j < species_num && i < cols.size(); j++)
		{
			orfs.back().align_evalues[j] = stod(cols[i]);
			i++;
		}
		orfs.back().aaseq = cols[i];
		i++;
		for (int j = 0; j < species_num && i < cols.size(); j++)
		{
			orfs.back().alt_aaseq[j] = cols[i];
			i++;
		}
		for (int j = 0; j < species_num && i < cols.size(); j++)
		{
			orfs.back().focal_align_aaseq[j] = cols[i];
			i++;
		}
		for (int j = 0; j < species_num && i < cols.size(); j++)
		{
			orfs.back().alt_align_aaseq[j] = cols[i];
			i++;
		}
		for (int k = 0; k < cols[i].size(); k++)
		{
			orfs.back().seq.push_back(nuc_map.at(cols[i][k]));
		}
		i++;
		for (int j = 0; j < species_num && i < cols.size(); j++)
		{
			for (int k = 0; k < cols[i].size(); k++)
			{
				orfs.back().alt_seqs[j].push_back(nuc_map.at(cols[i][k]));
			}
			i++;
		}
		orfs.back().repeat_count = stoi(cols[i++]);
		orfs.back().orf_class = cols[i++];

		for (int j = 0; j < species_num && i < cols.size(); j++)
		{
			for (int k = 0; k < cols[i].size(); k++)
			{
				orfs.back().focal_flank_align[j].push_back(nuc_map.at(cols[i][k]));
			}
			i++;
		}


		for (int j = 0; j < species_num && i < cols.size(); j++)
		{
			for (int k = 0; k < cols[i].size(); k++)
			{
				orfs.back().flank_align[j].push_back(nuc_map.at(cols[i][k]));
			}
			i++;
		}
		/*if (i >= cols.size())
		{
			continue;
		}
		orfs.back().start_codon_overlap = stoi(cols[i++]);
		orfs.back().stop_codon_overlap = stoi(cols[i++]);*/

		for (int j = 0; j < species_num; j++)
		{
			orfs.back().use_blast_align[j] = stoi(cols[i++]);
		}		
		for (int j = 0; j < species_num; j++)
		{
			orfs.back().homology_validated[j] = stoi(cols[i++]);
		}

	}
}

void read_orfs(vector<Orf>& orfs, string filename, int species_num)
{
	ifstream file(filename);
	string line;
	getline(file, line);
	while (getline(file, line))
	{
		int i = 1;
		vector<string> cols;
		split(line, ' ', cols);
		orfs.push_back(Orf(species_num));
		orfs.back().orf_id = cols[i++];
		orfs.back().contig = stoi(cols[i++]);
		orfs.back().start_pos = stoi(cols[i++]);
		orfs.back().end_pos = stoi(cols[i++]);
		orfs.back().strand = stoi(cols[i++]);
		orfs.back().overlaps_gene = cols[i++];
		orfs.back().overlap_gene_length = stoi(cols[i++]);
		orfs.back().overlap_gene_relative_frame = stoi(cols[i++]);
		orfs.back().is_gene = cols[i++];
		orfs.back().splice_gene = cols[i++];
		orfs.back().up_extend_gene = cols[i++];
		orfs.back().tifseqs_hits0 = stoi(cols[i++]);
		orfs.back().tifseqs_hits1 = stoi(cols[i++]);
		orfs.back().orf_class = cols[i++];
	}
}

void print_orf_overlap(const vector<Overlap>& overlaps, string filename)
{
	ofstream file(filename);
	file << "id overlap overlap_length overlap_relframe all_overlaps all_overlaps_length all_overlaps_relframe sense_ver sense_unchar sense_dub sense_pseudo sense_te sense_blocked anti_ver anti_unchar anti_dub anti_pseudo anti_te anti_blocked tRNA ltr sno ars sn ncRNA";
	for (int i = 0; i < overlaps.size(); i++)
	{
		if (overlaps[i].all_overlaps.size() > 0)
		{
			file << "\n" << i << " " << overlaps[i].greatest_overlap << " " << overlaps[i].greatest_overlap_length << " " << overlaps[i].greatest_overlap_rel_frame << " ";
			for (int j = 0; j < overlaps[i].all_overlaps.size(); j++)
			{
				file << overlaps[i].all_overlaps[j] << ",";
			}
			file << " ";
			for (int j = 0; j < overlaps[i].all_overlap_lengths.size(); j++)
			{
				file << overlaps[i].all_overlap_lengths[j] << ",";
			}
			file << " ";
			for (int j = 0; j < overlaps[i].all_overlap_rel_frames.size(); j++)
			{
				file << overlaps[i].all_overlap_rel_frames[j] << ",";
			}
			file << " " << overlaps[i].sense_verified_overlap << " " << overlaps[i].sense_uncharacterized_overlap << " " << overlaps[i].sense_dubious_overlap << " " << overlaps[i].sense_pseudogene_overlap << " " << overlaps[i].sense_te_overlap << " " << overlaps[i].sense_blocked_overlap;
			file << " " << overlaps[i].anti_verified_overlap << " " << overlaps[i].anti_uncharacterized_overlap << " " << overlaps[i].anti_dubious_overlap << " " << overlaps[i].anti_pseudogene_overlap << " " << overlaps[i].anti_te_overlap << " " << overlaps[i].anti_blocked_overlap;
			file << " " << overlaps[i].trna_overlap << " " << overlaps[i].ltr_overlap << " " << overlaps[i].sno_overlap << " " << overlaps[i].ars_overlap << " " << overlaps[i].sn_overlap << " " << overlaps[i].ncrna_overlap;
		}
		else
		{
			file << "\n" << i << " X 0 0 X X X 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
		}
	}
}

void print_shared_transcripts(const vector<SharedTranscript>& shared, string filename)
{
	ofstream file(filename);
	file << "id orf_id relation sgd_id sgd_class";
	for (int i = 0; i < shared.size(); i++)
	{
		file << "\n" << i << " " << shared[i].orf_id << " " << shared[i].relation << " " << shared[i].shared_sgd << " " << shared[i].shared_sgd_class;
	}
}

void read_riboseq_studies(vector<Experiment>& experiments, string filename)
{
	map<char, int> nucmap = { { 'A',0 },{ 'a',0 },{ 'C',1 },{ 'c',1 },{ 'G',2 },{ 'g',2 },{ 'T',3 },{ 't',3 },{ 'N',5 },{ 'n',5 },{'-',4} };
	ifstream file(filename);
	string line;
	getline(file, line);
	int index = 1;
	while (getline(file, line))
	{
		vector<string> columns;
		split(line, '\t', columns);
		experiments.push_back(Experiment());
		experiments.back().srr = columns[2];
		experiments.back().srp = columns[3];
		experiments.back().preprocessed = stoi(columns[4]);
		if (!experiments.back().preprocessed)
		{
			for (int i = 0; i < columns[5].size(); i++)
			{
				experiments.back().linker.push_back(nucmap.at(columns[5][i]));
			}
		}
		if (columns[6] == "Yes")
		{
			experiments.back().chx = 1;
		}
		else if (columns[6] == "No")
		{
			experiments.back().chx = 0;
		}
		if (columns[7] == "Yes")
		{
			experiments.back().ypd = 1;
		}
		else if (columns[7] == "No")
		{
			experiments.back().ypd = 0;
		}
		if (columns[11] == "Yes")
		{
			experiments.back().ypd = 4;
		}
	}
}

void read_ribseqs(vector<Ribseq>& rseqs, Experiment& experiment, int limit, int paired, int end_cutoff)
{
	map<char, int> nucmap = { { 'A',0 },{ 'a',0 },{ 'C',1 },{ 'c',1 },{ 'G',2 },{ 'g',2 },{ 'T',3 },{ 't',3 },{ 'N',4 },{ 'n',4 } };
	cout << "\nreading file: " << "SRR/SRR" + experiment.srr + ".fastq";
	ifstream file;//("SRR/"+experiment.srr + ".fastq"); //SRR2829333.fastq
	if (!paired)
	{
		file.open("SRR/" + experiment.srr + ".fastq");
	}
	else if (paired == 1)
	{
		file.open("SRR/" + experiment.srr + "_" + to_string(paired) + ".fastq");
	}
	string line;
	int i = 0;
	while (getline(file, line))
	{
		rseqs.push_back(Ribseq());
		vector<string> columns;
		split(line, ' ', columns);
		rseqs.back().id = columns[0];
		getline(file, line);
		for (int i = 0; i < line.size(); i++)
		{
			rseqs.back().seq.push_back(nucmap.at(line[i]));
		}
		rseqs.back().rc_seq = rseqs.back().seq;
		reverse_complement(rseqs.back().rc_seq);
		getline(file, line);
		getline(file, line);
		for (int i = 0; i < line.size(); i++)
		{
			rseqs.back().phred.push_back(line[i] - '!');
		}
		if (limit != -1 && rseqs.size() == limit)
		{
			break;
		}
	}
}

unsigned long int vector_to_hashval(const vector<int>& vec)
{
	unsigned long int hash_value = 0;
	for (int i = 0; i < vec.size(); i++)
	{
		hash_value += (unsigned long int)vec[i] * (unsigned long int)pow((unsigned long int)4, (unsigned long int)vec.size() - (unsigned long int)i - (unsigned long int)1);
	}
	return hash_value;
}

void filter_ribseqs(vector<Ribseq>& rseqs, const int min_size_class, const int max_size_class, const vector<int>& linker)
{
	vector<Ribseq> filtered_rseqs;
	int linker_length = 6;
	if (linker.size() < linker_length)
	{
		linker_length = linker.size();
	}
	vector<int> linker_used(linker.begin(), linker.begin() + linker_length);
	int kmer_length = linker_used.size();
	unsigned long int hash_linker = vector_to_hashval(linker_used);
	for (int i = 0; i < rseqs.size(); i++)
	{
		int current_length = 0;
		unsigned long int hashval = 0;
		bool qc_pass = true;
		int linker_pos = -1;
		bool rc = false;
		for (int j = 0; j < rseqs[i].seq.size(); j++)
		{
			if (rseqs[i].seq[j] < 4 && rseqs[i].phred[j] >= 20)
			{
				if (current_length < kmer_length)
				{
					hashval += (unsigned long int)rseqs[i].seq[j] * (unsigned long int)pow((unsigned long int)4, (unsigned long int)kmer_length - (unsigned long int)current_length - 1);
					current_length++;
				}
				else if (current_length == kmer_length)
				{
					hashval -= (unsigned long int)rseqs[i].seq[j - kmer_length] * (unsigned long int)pow(4, (unsigned long int)kmer_length - 1);
					hashval *= 4;
					hashval += (unsigned long int)rseqs[i].seq[j];
					if (hashval == hash_linker)
					{
						linker_pos = j + 1 - kmer_length;
						break;
					}
				}
			}
			else
			{
				qc_pass = false;
				break;
			}
		}

		if (qc_pass && linker_pos > -1)
		{
			vector<int> proc_seq;
			if (!rc)
			{
				proc_seq = vector<int>(rseqs[i].seq.begin() /*+ 1*/, rseqs[i].seq.begin() + linker_pos);
			}
			else
			{
				proc_seq = vector<int>(rseqs[i].rc_seq.begin() + 1, rseqs[i].rc_seq.begin() + linker_pos);
			}

			if ((proc_seq.size() >= min_size_class && proc_seq.size() <= max_size_class) || min_size_class == -1)
			{
				filtered_rseqs.push_back(rseqs[i]);
				filtered_rseqs.back().proc_seq = proc_seq;
				filtered_rseqs.back().hash = vector_to_hashval(proc_seq);
			}
		}
	}
	rseqs = filtered_rseqs;
}

void cut_ribseqs(vector<Ribseq>& rseqs, int end_cutoff, int start_cutoff)
{
	for (int i = 0; i < rseqs.size(); i++)
	{
		rseqs[i].proc_seq = vector<int>(rseqs[i].proc_seq.begin() + start_cutoff, rseqs[i].proc_seq.end() - end_cutoff);
		rseqs[i].hash = vector_to_hashval(rseqs[i].proc_seq);
	}
}

void filter_preprocessed_ribseqs(vector<Ribseq>& rseqs, const int min_size_class, const int max_size_class, int remove_nuc5)
{
	vector<Ribseq> filtered_rseqs;
	for (int i = 0; i < rseqs.size(); i++)
	{
		bool qc_pass = true;
		for (int j = 0; j < rseqs[i].seq.size(); j++)
		{
			if (rseqs[i].seq[j] < 4 && rseqs[i].phred[j] >= 20)
			{
			}
			else
			{
				qc_pass = false;
			}
		}
		if (qc_pass && ((rseqs[i].seq.size() >= min_size_class && rseqs[i].seq.size() <= max_size_class) || min_size_class == -1))
		{
			filtered_rseqs.push_back(rseqs[i]);
			filtered_rseqs.back().proc_seq = vector<int>(filtered_rseqs.back().seq.begin() + remove_nuc5, filtered_rseqs.back().seq.end());
			filtered_rseqs.back().hash = vector_to_hashval(filtered_rseqs.back().proc_seq);
		}
	}
	rseqs = filtered_rseqs;
}

void get_rseqs_by_length(vector<Ribseq>& rseqs0, const vector<Ribseq>& rseqs, const int length)
{
	for (int i = 0; i < rseqs.size(); i++)
	{
		if (rseqs[i].proc_seq.size() == length /*- 1*/)
		{
			rseqs0.push_back(rseqs[i]);
		}
	}
}

void collect_hashes(map<long unsigned int, vector<int> >& hashes, const vector<Ribseq>& rseqs)
{
	for (int i = 0; i < rseqs.size(); i++)
	{
		hashes[rseqs[i].hash].push_back(i);
	}
}

void map_ribseqs(vector<Ribseq>& rseqs, const map<long unsigned int, vector<int> >& hashes, const vector<Contig>& contigs, int kmer_length, int strand)
{
	for (int i = 0; i < contigs.size(); i++)
	{
		int current_length = 0;
		unsigned long int hashval = 0;
		for (int j = 0; j < contigs[i].seq.size(); j++)
		{
			if (contigs[i].seq[j] < 4)
			{
				if (current_length < kmer_length)
				{
					hashval += (unsigned long int)contigs[i].seq[j] * (unsigned long int)pow((unsigned long int)4, (unsigned long int)kmer_length - (unsigned long int)current_length - (unsigned long int)1);
					current_length++;
				}
				else if (current_length == kmer_length)
				{
					hashval -= (unsigned long int)contigs[i].seq[j - kmer_length] * (unsigned long int)pow((unsigned long int)4, (unsigned long int)kmer_length - (unsigned long int)1);
					hashval *= 4;
					hashval += (unsigned long int)contigs[i].seq[j];
					if (hashes.count(hashval))
					{
						for (int k = 0; k < hashes.at(hashval).size(); k++)
						{
							Pos pos;
							pos.contig = i;
							pos.coord = j - kmer_length;
							pos.strand = strand;
							if (strand == 1)
							{
								pos.coord = contigs[i].seq.size() - j - 1 + kmer_length;
							}
							rseqs.at(hashes.at(hashval)[k]).mappings.push_back(pos);
						}
					}
				}
			}
		}
	}
}

void filter_mapped_ribseqs(vector<Ribseq>& rseqs)
{
	vector<Ribseq> filtered_rseqs;
	for (int i = 0; i < rseqs.size(); i++)
	{
		if (rseqs[i].mappings.size() == 1)
		{
			filtered_rseqs.push_back(rseqs[i]);
		}
	}
	rseqs = filtered_rseqs;
}

void place_mapped_ribseqs(vector<vector<int>>& mapped_ribseqs, const vector<Ribseq>& rseqs, const vector<Contig>& contigs_focal)
{
	mapped_ribseqs = vector<vector<int>>(contigs_focal.size());
	for (int i = 0; i < contigs_focal.size(); i++)
	{
		mapped_ribseqs[i] = vector<int>(contigs_focal[i].seq.size());
	}
	for (int i = 0; i < rseqs.size(); i++)
	{
		if (rseqs[i].mappings.size() > 0 /*&& rseqs[i].mappings.back().strand == 0*/)
		{
			mapped_ribseqs[rseqs[i].mappings.back().contig][rseqs[i].mappings.back().coord]++;
		}
	}
}

void print_ribseqs(const vector<Ribseq>& rseqs, int rseq_length, string suffix)
{
	map<int, char> rev_nuc_map = { { 0,'A' },{ 1,'C' },{ 2,'G' },{ 3,'T' },{ 4,'N' },{ 5,'N' },{ 6,'N' } };

	ofstream file("riboseq_output_complete/ribseqs_" + suffix);
	file << "id full_length length hits contig strand pos hash seq";
	for (int i = 0; i < rseqs.size(); i++)
	{
		file << "\n" << i << " " << rseqs[i].seq.size() << " " << rseqs[i].proc_seq.size() << " " << rseqs[i].mappings.size();// << " " << rseqs[i].mappings;
		if (rseqs[i].mappings.size() > 0)
		{
			file << " " << rseqs[i].mappings.back().contig << " " << rseqs[i].mappings.back().strand << " " << rseqs[i].mappings.back().coord;
		}
		else
		{
			file << " -1 -1 -1";
		}
		file << " " << rseqs[i].hash;
		file << " ";
		for (int j = 0; j < rseqs[i].proc_seq.size(); j++)
		{
			file << rev_nuc_map.at(rseqs[i].proc_seq[j]);
		}
	}
}

void map_ribseqs_to_orfs(vector<Orf>& orfs, const vector<vector<int>>& mapped_ribseqs, int direction)
{
	for (int i = 0; i < orfs.size(); i++)
	{
		if (orfs[i].strand == 0 && direction == 0)
		{
			int frame = 0;
			for (int j = 0; j < orfs[i].seq.size(); j++)
			{
				orfs[i].frame_hits[frame] += (mapped_ribseqs[orfs[i].contig][orfs[i].start_pos + j] > 0);
				orfs[i].sum_hits[frame] += mapped_ribseqs[orfs[i].contig][orfs[i].start_pos + j];
				frame++;
				if (frame == 3)
				{
					frame = 0;
				}
			}
		}
		else if (orfs[i].strand == 1 && direction == 1)
		{
			int frame = 0;
			for (int j = 0; j < orfs[i].seq.size(); j++)
			{
				orfs[i].frame_hits[frame] += (mapped_ribseqs[orfs[i].contig][orfs[i].end_pos - j] > 0);
				orfs[i].sum_hits[frame] += mapped_ribseqs[orfs[i].contig][orfs[i].end_pos - j];
				frame++;
				if (frame == 3)
				{
					frame = 0;
				}
			}
		}
	}
}

int determine_p_site(const vector<Orf>& orfs, const vector<vector<int>> mapped_ribseqs)
{
	vector<int> rel_hits(40);
	for (int i = 0; i < orfs.size(); i++)
	{
		if (orfs[i].strand == 0)
		{
			for (int j = 0; j < rel_hits.size(); j++)
			{
				if (orfs[i].start_pos - 20 + j >= 0 && orfs[i].start_pos - 20 + j < mapped_ribseqs[orfs[i].contig].size())
				{
					rel_hits[j] += mapped_ribseqs[orfs[i].contig][orfs[i].start_pos - 20 + j];
				}
			}
		}
	}
	int sum_rel_hits = 0;
	vector<int> frame_hits(3);
	for (int i = 0; i < rel_hits.size(); i++)
	{
		sum_rel_hits += rel_hits[i];
		frame_hits[i % 3] += rel_hits[i];
	}
	int best_frame = 2;
	if (frame_hits[0] > frame_hits[1] && frame_hits[0] > frame_hits[2])
	{
		best_frame = 0;
	}
	else if (frame_hits[1] > frame_hits[2])
	{
		best_frame = 1;
	}
	int p = 0;
	for (int i = 0; i < rel_hits.size(); i++)
	{
		if ((double)rel_hits[i] / (double)sum_rel_hits > .05 && i % 3 == best_frame)
		{
			p = i;
			break;
		}
	}
	return p - 20;
}

void replace_mapped_ribseqs(vector<vector<int>>& mapped_ribseqs, const vector<Ribseq>& rseqs, const vector<Contig>& contigs_focal, int p_site, int direction)
{
	mapped_ribseqs = vector<vector<int>>(contigs_focal.size());
	for (int i = 0; i < contigs_focal.size(); i++)
	{
		mapped_ribseqs[i] = vector<int>(contigs_focal[i].seq.size());
	}
	for (int i = 0; i < rseqs.size(); i++)
	{
		if (rseqs[i].mappings.size() > 0 && rseqs[i].mappings.back().strand == direction)
		{
			if (direction == 0)
			{
				mapped_ribseqs[rseqs[i].mappings.back().contig][rseqs[i].mappings.back().coord - p_site]++;
			}
			else
			{
				mapped_ribseqs[rseqs[i].mappings.back().contig][rseqs[i].mappings.back().coord + p_site]++;
			}
		}
	}
}

void count_frame_hits(vector<int>& frame_hits, const vector<Orf>& orfs, int strand)
{
	frame_hits = vector<int>(3);
	for (int i = 0; i < orfs.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (orfs[i].strand == strand || strand == -1)
			{
				frame_hits[j] += orfs[i].frame_hits[j];
			}
		}
	}
}

void print_mapped_ribseqs(const vector<vector<int>>& mapped_ribseqs, vector<Contig>& contigs_focal, int rseq_length, string suffix)
{
	ofstream file("riboseq_output_complete/mapped_ribseqs_" + suffix);
	file << "contig pos nuc count";
	for (int i = 0; i < mapped_ribseqs.size(); i++)
	{
		for (int j = 0; j < mapped_ribseqs[i].size(); j++)
		{
			file << "\n" << i << " " << j << " " << contigs_focal[i].seq[j] << " " << mapped_ribseqs[i][j];
		}
	}
}



void map_riboseq_reads(int start_map, int end_map)
{
	//maps riboseq reads to genome
	vector<Orf> orfs;
	read_annotated_genes(orfs, "input_files/orf_coding.fasta", true);//read sgd ORF annotations
	cout << "\nannotated genes read: " << orfs.size();
	sort(orfs.begin(), orfs.end());
	vector<Contig> contigs_focal;
	read_fasta(contigs_focal, "input_files/S288C_reference_sequence_R64-2-1_20150113.fsa");//read s cerevisiae genome
	vector<Contig> rc_contigs_focal = contigs_focal;
	rev_comp_contigs(rc_contigs_focal);//create reverse complemented contigs of the yeast genome so we can analyze the same way in either direction
	cout << "\nread genome";

	vector<Experiment> experiments;
	vector<int> rseq_lengths = { 25,26,27,28,29,30,31,32,33,34,35 };//lengths of riboseq reads to analyze
	read_riboseq_studies(experiments, "input_files/riboseq_experiments_dataset.txt");//this file lists information on the fastq files to read
	cout << "\nread riboseq studies: " << experiments.size();

	for (int i = start_map; i <= end_map; i++) //sets the range of experiments to process. to analyze all, set start_map=0, end_map=experiments.size()
	{
		ofstream file("riboseq_output_complete/riboseqs_dataset_" + experiments[i].srr); //files for output data and QC
		ofstream info_file("riboseq_output_complete/riboseqs_infofile_" + experiments[i].srr);
		file << "filename study id length frame0 frame1 frame2 p_site";
		info_file << "total_reads filtered_reads";
		for (int k = 0; k < rseq_lengths.size(); k++)
		{
			info_file << " length" << rseq_lengths[k] << " mapped_length" << rseq_lengths[k] << " unique_length" << rseq_lengths[k];
		}
		info_file << "\n";
		cout << "\nprocessing experiment " << experiments[i].srr << " preprocessed: " << experiments[i].preprocessed << " adaptor: " << experiments[i].linker.size();
		int total_mapped = 0;
		vector<Ribseq> rseqs;
		if (experiments[i].paired == 0)
		{
			read_ribseqs(rseqs, experiments[i], -1, 0, experiments[i].end_cutoff);
		}
		else
		{
			read_ribseqs(rseqs, experiments[i], -1, 1, experiments[i].end_cutoff);
			read_ribseqs(rseqs, experiments[i], -1, 2, experiments[i].end_cutoff);
		}
		cout << "\nribseqs read: " << rseqs.size();
		info_file << rseqs.size() << " ";
		if (!experiments[i].preprocessed)//if the reads are not preprocessed, we need to find and trim the adapters
		{
			filter_ribseqs(rseqs, rseq_lengths.front(), rseq_lengths.back(), experiments[i].linker);
		}
		else
		{
			filter_preprocessed_ribseqs(rseqs, rseq_lengths.front(), rseq_lengths.back(), 1);
		}
		cut_ribseqs(rseqs, experiments[i].end_cutoff, experiments[i].start_cutoff);//trim adapters
		cout << "\nribseqs filtered: " << rseqs.size();
		info_file << rseqs.size() << " ";
		for (int k = 0; k < rseq_lengths.size(); k++)//for each experiment, each ribo-seq read length is considered individually 
		{
			file << "\n" << experiments[i].srr << " " << experiments[i].srp << " " << rseq_lengths[k] << " ";
			cout << "\n" << experiments[i].srr << " " << experiments[i].srp << " " << rseq_lengths[k] << " ";

			vector<Ribseq> rseqs0;
			get_rseqs_by_length(rseqs0, rseqs, rseq_lengths[k]);
			info_file << rseqs0.size() << " ";
			cout << "\nreads at length " << rseq_lengths[k] << ": " << rseqs0.size();
			map<long unsigned int, vector<int>> hashes;
			collect_hashes(hashes, rseqs0);
			map_ribseqs(rseqs0, hashes, contigs_focal, rseq_lengths[k], 0);//reads are mapped to the genome on both strands
			map_ribseqs(rseqs0, hashes, rc_contigs_focal, rseq_lengths[k], 1);
			info_file << rseqs0.size() << " ";
			string suffix0 = experiments[i].srr + "_" + to_string(rseq_lengths[k]) + "_all";
			filter_mapped_ribseqs(rseqs0);
			info_file << rseqs0.size() << " ";

			total_mapped += rseqs0.size();
			cout << "\nribseqs uniquely mapped: " << rseqs0.size() << " " << total_mapped;
			if (rseqs0.size() > 10000)//we only consider the read length further if there are at least 10000 reads at that length
			{
				vector<vector<int>> mapped_ribseqs;
				place_mapped_ribseqs(mapped_ribseqs, rseqs0, contigs_focal);
				string suffix = experiments[i].srr + "_" + to_string(rseq_lengths[k]);
				print_ribseqs(rseqs0, rseq_lengths[k], suffix);

				vector<Orf> orfs0 = orfs;
				map_ribseqs_to_orfs(orfs0, mapped_ribseqs, 0);
				int p_site = determine_p_site(orfs0, mapped_ribseqs);//determine the p-site associated with the read length based on the read pattern of reads near the start codons of annotated genes

				vector<vector<int>> mapped_forward_ribseqs;
				vector<vector<int>> mapped_reverse_ribseqs;

				replace_mapped_ribseqs(mapped_forward_ribseqs, rseqs0, contigs_focal, p_site, 0);//remap reads to inferred p-site position (rather than start position of the read). note that this is a simple uniform shift in position from the original mapping. in this way all read lengths will be comparable and reads should start at the actual start codon
				replace_mapped_ribseqs(mapped_reverse_ribseqs, rseqs0, contigs_focal, p_site, 1);

				//output data on how often reads map to firt, second, third position in codon. This is used to assess whether this read length in this experiment exhibits triplet periodicity.
				int test = determine_p_site(orfs0, mapped_forward_ribseqs);
				vector<int> frame_hits(3);
				count_frame_hits(frame_hits, orfs0, -1);
				print_mapped_ribseqs(mapped_forward_ribseqs, contigs_focal, rseq_lengths[k], suffix + "_f");
				print_mapped_ribseqs(mapped_reverse_ribseqs, contigs_focal, rseq_lengths[k], suffix + "_r");

				file << frame_hits[0] << " " << frame_hits[1] << " " << frame_hits[2] << " " << p_site;
			}
		}
	}
}

void print_orf_neighbors(const vector<Orf>& orfs, string filename)
{
	ofstream file(filename);
	file << "id up_neighbor up_neighbor_dist down_neighbor down_neighbor_dist start_codon_overlap stop_codon_overlap";
	for (int i = 0; i < orfs.size(); i++)
	{
		file << "\n" << i << " " << orfs[i].up_gene_sys << " " << orfs[i].up_gene_dist << " " << orfs[i].down_gene_sys << " " << orfs[i].down_gene_dist<<" "<<orfs[i].start_codon_overlap<<" "<<orfs[i].stop_codon_overlap;
	}
}

void get_orf_neighbors(vector<Orf>& orfs, const vector<SGD>& genes)
{
	int start_orf = 0;
	for (int i = 0; i < genes.size() - 1; i++)
	{
		cout << "\ngene " << i;
		while (orfs[start_orf].contig < genes[i].chr)
		{
			start_orf++;
		}
		for (int j = start_orf; j < orfs.size(); j++)
		{
			if (orfs[j].start_pos > genes[i + 1].end)
			{
				break;
			}
			if (orfs[j].contig == genes[i].chr && genes[i].chr == genes[i + 1].chr)
			{
				if (orfs[j].start_pos == genes[i].start)
				{
					orfs[j].up_gene_sys = genes[i - 1].Name;
					orfs[j].down_gene_sys = genes[i + 1].Name;
					orfs[j].up_gene_dist = orfs[j].start_pos - genes[i - 1].end;
					orfs[j].down_gene_dist = genes[i + 1].start - orfs[j].end_pos;
				}
				if (orfs[j].start_pos > genes[i].start&& orfs[j].end_pos < genes[i + 1].end)
				{
					orfs[j].up_gene_sys = genes[i].Name;
					orfs[j].down_gene_sys = genes[i + 1].Name;
					orfs[j].up_gene_dist = orfs[j].start_pos - genes[i].end;
					orfs[j].down_gene_dist = genes[i+1].start - orfs[j].end_pos;
				}
				if(orfs[j].strand==1)
				{
					string dummy=orfs[j].up_gene_sys;
					orfs[j].up_gene_sys=orfs[j].down_gene_sys;
					orfs[j].down_gene_sys=dummy;
					int dummy_dist=orfs[j].up_gene_dist;
					orfs[j].up_gene_dist=orfs[j].down_gene_dist;
					orfs[j].down_gene_dist=dummy_dist;
				}
			}
		}
	}
}

void get_all_orfs()
{
	//download SGD annotations and acquire list of ORFs from scanning the genome 
	
	string suffix = "";
	vector<Contig> contigs; //since the genome we are using is assembled, each contig is an S. cerevisiae chromosome, 0-indexed
	vector<SGD> sgds; //SGD annotations

	read_sgd_annotation(sgds, "input_files/saccharomyces_cerevisiae.gff"); //annotation file downloaded from SGD
	cout << "\nread SGD annotations: " << sgds.size();

	vector<SGD> sgds_organized;
	organize_sgds(sgds_organized, sgds); //reorganize input SGD
	sort(sgds_organized.begin(), sgds_organized.end());
	cout << "\nSGD annotations organized: " << sgds_organized.size();

	vector<Orf> genes; 
	read_annotated_genes(genes, "input_files/orf_coding_all.fasta", false); //all annotated coding ORFs downloaded from SGD
	cout << "\nannotated genes read: " << genes.size();
	read_annotated_genes(genes, "input_files/rna_coding.fasta", true); //RNA genes
	read_annotated_genes(genes, "input_files/other_features_genomic.fasta", false);//other annotations on SGD
	sort(genes.begin(), genes.end());

	vector<Orf> orfs;
	read_fasta(contigs, "input_files/S288C_reference_sequence_R64-2-1_20150113.fsa"); //read genome downloaded from SGD
	get_orfs(orfs, contigs, 1, 7, false); //identify all single-exon ORFs in genome, start-to-stop (ATG starts) 

	remove_duplicate_orfs_by_contig(orfs, contigs); //any set of ORFs that share a stop codon are considered the same ORF. Remove all but the longest
	cout << "\norfs remaining after remove duplicates: " << orfs.size();

	//determine which annotations each ORF overlaps 
	
	cout << "\ncheck orf overlap";
	vector<Overlap> overlaps;

	sort(orfs.begin(), orfs.end());
	check_orf_overlap(orfs, sgds_organized, overlaps);
	get_orf_neighbors(orfs, sgds_organized);//determine closest annotated genes for each ORF 

	//determine which known transcript each ORF is located on, if any. Based on tifseq data from Pelechano et al. 2014. Also determine if ORFs share a transcript with any annotatd genes
	vector<Tifseq> tseqs;
	vector<SharedTranscript> shared;
	read_tifseqs(tseqs, "input_files/tsedAnno_V2.txt"); //taken from Pelechano et al. 2014 supplementary data
	cout << "\ntifseqs read: " << tseqs.size();
	sort(tseqs.begin(), tseqs.end());

	get_tifseq_coverage(orfs, tseqs, false);
	get_tifseq_coverage(sgds_organized, tseqs, true);
	get_tifseq_context(shared, orfs, sgds_organized, tseqs);

	sort_orfs_by_strand(orfs, overlaps);

	print_orfs_aa_fasta("orfs_comp" + suffix + ".fasta", orfs, true); //print amino acid sequence of each ORF in fasta format
	print_orfs_aa_fasta_scrambled_and_real("orfs_comp_scrambled_and_real" + suffix + ".fasta", orfs, true); //print amino acid sequence of each ORF, together with scrambled sequences of each ORF to use as a negative control

	print_orfs_nuc_fasta("orfs_comp_nuc" + suffix + ".fasta", orfs);//print nucleotide sequence of each ORF

	//output data tables for identified ORFs 
	print_orfs(orfs, "_comp" + suffix);
	print_orf_overlap(overlaps, "orf_overlap" + suffix);
	print_shared_transcripts(shared, "shared_transcripts" + suffix);
	print_orf_neighbors(orfs, "orf_neighbors");
}

void combo_ribosq_qc(vector<Experiment> experiments)
{
	ofstream file("riboseqs_dataset_combo");
	file << "filename study id length frame0 frame1 frame2 p_site";
	for (int i = 0; i < experiments.size(); i++)
	{
		ifstream infile("riboseq_output_complete/riboseqs_dataset_" + experiments[i].srr);
		string line;
		getline(infile, line);
		while (getline(infile, line))
		{
			file << "\n" << line;
		}
	}
}

void get_riboseq_qc(vector<Experiment>& experiments, string filename, bool force_pass)
{
	map<string, int> exp_map;
	for (int i = 0; i < experiments.size(); i++)
	{
		exp_map[experiments[i].srr] = i;
	}

	ifstream file(filename);
	string line;
	getline(file, line);
	while (getline(file, line))
	{
		vector<string> column_data;
		split(line, ' ', column_data);
		string srr = column_data[0];
		int read_length = stoi(column_data[2]);
		if (column_data.size() > 4)
		{
			vector<int> frames = { stoi(column_data[3]),stoi(column_data[4]) ,stoi(column_data[5]) };
			for (int i = 0; i < frames.size(); i++)
			{
				vector<int> other_frames;
				for (int j = 0; j < frames.size(); j++)
				{
					if (j != i)
					{
						other_frames.push_back(j);
					}
				}
				bool allgood = true;
				for (int j = 0; j < other_frames.size(); j++)
				{
					if (frames[i] < frames[other_frames[j]] * 2)
					{
						allgood = false;
					}
				}
				if (allgood && frames[i] > 10000 || force_pass)
				{
					experiments[exp_map.at(srr)].frames[read_length] = i;
				}
			}
		}
	}
}

void read_mapped_ribseqs(vector<vector<int>>& mapped_ribseqs, string filename)
{
	ifstream file(filename);
	string line;
	getline(file, line);
	int index = 1;
	while (getline(file, line))
	{
		vector<string> columns;
		split(line, ' ', columns);
		int contig = stoi(columns[0]);
		if (contig == mapped_ribseqs.size())
		{
			mapped_ribseqs.push_back(vector<int>());
		}
		mapped_ribseqs.back().push_back(stoi(columns[3]));
	}
}

void print_mapped_ribseqs(const vector<vector<vector<int>>>& mapped_ribseqs, vector<Contig>& contigs_focal, int rseq_length, string suffix, vector<string> filenames)
{
	ofstream file("riboseq_output_complete/mapped_ribseqs_" + suffix);
	file << "contig pos nuc";
	for (int k = 0; k < mapped_ribseqs.size(); k++)
	{
		file << " " << filenames[k];
	}
	for (int i = 0; i < contigs_focal.size(); i++)
	{
		for (int j = 0; j < contigs_focal[i].seq.size(); j++)
		{
			file << "\n" << i << " " << j << " " << contigs_focal[i].seq[j];// << " " << mapped_ribseqs[i][j];
			for (int k = 0; k < mapped_ribseqs.size(); k++)
			{
				file << " " << mapped_ribseqs[k][i][j];
			}
		}
	}
}

void combine_riboseq_reads(int start_map,int end_map)	
{//assembles and combines mapped riboseq reads from a collection of studies
	vector<Contig> contigs_focal;
	read_fasta(contigs_focal, "input_files/S288C_reference_sequence_R64-2-1_20150113.fsa");

	vector<Experiment> experiments;
	vector<int> rseq_lengths = { 25,26,27,28,29,30,31,32,33,34,35 };
	read_riboseq_studies(experiments, "input_files/riboseq_experiments_dataset.txt");
	cout << "\nread riboseq studies: " << experiments.size() << " " << experiments.back().srr << " " << experiments.front().srr;
	combo_ribosq_qc(experiments);//combines the individual riboseq qc files into a larger file
	get_riboseq_qc(experiments, "riboseqs_dataset_combo", false);//read qc data for riboseq studies
	cout << "\nget riboseq qc";
	vector<string> filenames;
	vector<int> frames;
	vector<vector<vector<int>>> mapped_all_forward;
	vector<vector<vector<int>>> mapped_all_reverse;
	if (start_map == -1)
	{
		start_map = 0;
		end_map = experiments.size() - 1;
	}
	//reads mapped riboseq reads from individual files corresponding to individual experiments and combines into a single file
	for (int i = start_map; i <= end_map; i++)
	{
		for (int k = 0; k < rseq_lengths.size(); k++)
		{
			if (experiments[i].frames.count(rseq_lengths[k]))
			{
				int frame = experiments[i].frames.at(rseq_lengths[k]);
				vector<vector<int>> mapped_ribseqs_forward;
				vector<vector<int>> mapped_ribseqs_reverse;

				string filename = "riboseq_output_complete/mapped_ribseqs_" + experiments[i].srr + "_" + to_string(rseq_lengths[k]);
				cout << "\nprocessing " << filename;
				read_mapped_ribseqs(mapped_ribseqs_forward, filename + "_f");
				read_mapped_ribseqs(mapped_ribseqs_reverse, filename + "_r");

				mapped_all_forward.push_back(mapped_ribseqs_forward);
				mapped_all_reverse.push_back(mapped_ribseqs_reverse);

				filenames.push_back(experiments[i].srr + "_" + to_string(rseq_lengths[k]));
			}
		}
	}
	print_mapped_ribseqs(mapped_all_forward, contigs_focal, 0, "combo_all_f", filenames);
	print_mapped_ribseqs(mapped_all_reverse, contigs_focal, 0, "combo_all_r", filenames);
}
void read_riboseq_reads(vector<Ribopos>& ribopos, string filename, int strand, int readlength, int specific_study, vector<string>& studies)
{
	studies = vector<string>();
	ifstream file(filename);
	vector<string> cols;
	string line;
	getline(file, line);
	split(line, ' ', cols);
	vector<int> length_matches;
	if (specific_study == -1)
	{
		for (int i = 3; i < cols.size(); i++)
		{
			if (readlength == -1 || stoi(cols[i].substr(cols[i].size() - 2, 2)) == readlength)
			{
				length_matches.push_back(i);
				int pos = cols[i].find("_");
				studies.push_back(cols[i].substr(0, pos));
			}
		}
	}
	else
	{
		length_matches.push_back(3 + specific_study);
	}
	while (getline(file, line))
	{
		cols = vector<string>();
		split(line, ' ', cols);
		ribopos.push_back(Ribopos());
		ribopos.back().chr = stoi(cols[0]);
		ribopos.back().dir = strand;
		ribopos.back().pos = stoi(cols[1]);
		ribopos.back().nuc = stoi(cols[2]);
		ribopos.back().readlength = readlength;
		for (int i = 0; i < length_matches.size(); i++)
		{
			ribopos.back().readsum += stoi(cols[length_matches[i]]);
			ribopos.back().reads.push_back(stoi(cols[length_matches[i]]));
		}
	}
}

void map_riboseq_to_orfs(vector<Orf>& orfs, const vector<Ribopos>& ribopos, int strand, int specific_studies, set<int> studies_included)
{
	vector<map<int, vector<int>>> orfmap(17);
	vector<map<int, vector<int>>> orfmap_pos(17);
	for (int i = 0; i < orfs.size(); i++)
	{
		if (orfs[i].strand == strand)
		{
			orfs[i].riboseq_reads = vector<int>(1 + orfs[i].end_pos - orfs[i].start_pos);
			for (int j = orfs[i].start_pos; j <= orfs[i].end_pos; j++)
			{
				orfmap[orfs[i].contig][j].push_back(i);
				orfmap_pos[orfs[i].contig][j].push_back(j - orfs[i].start_pos);
			}
		}
	}
	for (int i = 0; i < ribopos.size(); i++)
	{
		if (orfmap[ribopos[i].chr].count(ribopos[i].pos))
		{
			for (int j = 0; j < orfmap[ribopos[i].chr][ribopos[i].pos].size(); j++)
			{
				if (specific_studies == -1)
				{
					orfs[orfmap[ribopos[i].chr][ribopos[i].pos][j]].riboseq_reads[orfmap_pos[ribopos[i].chr][ribopos[i].pos][j]] += ribopos[i].readsum;
					orfs[orfmap[ribopos[i].chr][ribopos[i].pos][j]].sum_mapped_reads += ribopos[i].readsum;
				}
				else
				{
					for (int k = 0; k < ribopos[i].reads.size(); k++)
					{
						if (studies_included.count(k))
						{
							orfs[orfmap[ribopos[i].chr][ribopos[i].pos][j]].riboseq_reads[orfmap_pos[ribopos[i].chr][ribopos[i].pos][j]] += ribopos[i].reads[k];
							orfs[orfmap[ribopos[i].chr][ribopos[i].pos][j]].sum_mapped_reads += ribopos[i].reads[k];
						}
					}
				}
			}
		}
	}
	if (strand == 1)
	{
		for (int i = 0; i < orfs.size(); i++)
		{
			if (orfs[i].strand == 1)
			{
				reverse(orfs[i].riboseq_reads.begin(), orfs[i].riboseq_reads.end());
			}
		}
	}
}

void scramble_riboseq_reads(vector<Orf>& orfs)
{
	for (int i = 0; i < orfs.size(); i++)
	{
		orfs[i].scrambled_reads = orfs[i].riboseq_reads;
		random_shuffle(orfs[i].scrambled_reads.begin(), orfs[i].scrambled_reads.end());
	}
}

double log_nChoosek(unsigned n, unsigned k)
{
	if (k == 0) return 0;

	double result = log(n);
	for (int i = 2; i <= k; ++i) {
		result += log(n - i + 1);
		result -= log(i);
	}
	return result;
}

double binom_test(int k, int n, double p)
{
	double pval = 0;
	for (int i = k; k <= n; k++)
	{
		double log_single_prob = log_nChoosek(n, k) + k * log(p) + (n - k) * log(1.0 - p);
		pval += exp(log_single_prob);
	}
	return pval;
}

void calc_riboseq_score(vector<riboseq_scorer>& scorer, const vector<Orf>& orfs)
{
	scorer = vector<riboseq_scorer>(orfs.size());
	for (int i = 0; i < orfs.size(); i++)
	{
		for (int j = 0; j < orfs[i].riboseq_reads.size(); j += 3)
		{
			int reads0 = orfs[i].riboseq_reads[j];
			int reads1 = orfs[i].riboseq_reads[j + 1];
			int reads2 = orfs[i].riboseq_reads[j + 2];
			if (reads0 > reads1&& reads0 > reads2)
			{
				scorer[i].hits[0]++;
			}
			else if (reads1 > reads0&& reads1 > reads2)
			{
				scorer[i].hits[1]++;
			}
			else if (reads2 > reads1&& reads2 > reads0)
			{
				scorer[i].hits[2]++;
			}
			scorer[i].frame_hits[0] += reads0;
			scorer[i].frame_hits[1] += reads1;
			scorer[i].frame_hits[2] += reads2;
		}
		for (int j = 0; j < orfs[i].scrambled_reads.size(); j += 3)
		{
			int reads0 = orfs[i].scrambled_reads[j];
			int reads1 = orfs[i].scrambled_reads[j + 1];
			int reads2 = orfs[i].scrambled_reads[j + 2];
			if (reads0 > reads1&& reads0 > reads2)
			{
				scorer[i].scrambled_hits[0]++;
			}
			else if (reads1 > reads0&& reads1 > reads2)
			{
				scorer[i].scrambled_hits[1]++;
			}
			else if (reads2 > reads1&& reads2 > reads0)
			{
				scorer[i].scrambled_hits[2]++;
			}
		}
		scorer[i].pval = binom_test(scorer[i].hits[0], scorer[i].hits[0] + scorer[i].hits[1] + scorer[i].hits[2], 1.0 / 3.0);
		scorer[i].scrambled_pval = binom_test(scorer[i].scrambled_hits[0], scorer[i].scrambled_hits[0] + scorer[i].scrambled_hits[1] + scorer[i].scrambled_hits[2], 1.0 / 3.0);
		if (false)//scorer[i].pval<.0001)
		{
			cout << "\norf: " << i << " " << scorer[i].pval << " " << scorer[i].hits[0] << " " << scorer[i].hits[1] << " " << scorer[i].hits[2];
			getchar();
		}
	}
}

void print_riboseqs_by_orf(const vector<riboseq_scorer>& scorer, const vector<Orf>& orfs, string suffix)
{
	ofstream file("riboseq_orfs" + suffix);
	file << "id reads scrambled hits0 hits1 hits2 scram0 scram1 scram2 frame0 frame1 frame2 pval scram_pval best_codon best_start best_pval second_codon second_start second_pval";
	for (int i = 0; i < orfs.size(); i++)
	{
		file << "\n" << i << " ";
		for (int j = 0; j < orfs[i].riboseq_reads.size(); j++)
		{
			file << orfs[i].riboseq_reads[j] << ",";
		}
		file << " ";
		for (int j = 0; j < orfs[i].scrambled_reads.size(); j++)
		{
			file << orfs[i].scrambled_reads[j] << ",";
		}
		file << " " << scorer[i].hits[0] << " " << scorer[i].hits[1] << " " << scorer[i].hits[2] << " " << scorer[i].scrambled_hits[0] << " " << scorer[i].scrambled_hits[1] << " " << scorer[i].scrambled_hits[2] << " " << scorer[i].frame_hits[0] << " " << scorer[i].frame_hits[1] << " " << scorer[i].frame_hits[2] << " " << scorer[i].pval << " " << scorer[i].scrambled_pval;
		file << " " << scorer[i].best_start_codon << " " << scorer[i].best_start_pos << " " << scorer[i].best_start_pval;
		file << " " << scorer[i].second_best_start_codon << " " << scorer[i].second_best_start_pos << " " << scorer[i].second_best_start_pval;
	}
}

void assemble_ribopos(vector<Ribopos>& combined_ribopos, const vector<Ribopos>& ribopos, const set<int>& included_columns)
{
	combined_ribopos = ribopos;
	for (int i = 0; i < ribopos.size(); i++)
	{
		combined_ribopos[i].reads[0] = 0;
		for (int j = 0; j < ribopos[i].reads.size(); j++)
		{
			if (included_columns.count(j))
			{
				combined_ribopos[i].reads[0] += ribopos[i].reads[j];
			}
		}
	}
}

int count_reads(const vector<Ribopos>& ribopos)
{
	int count = 0;
	for (int i = 0; i < ribopos.size(); i++)
	{
		count += ribopos[i].reads[0];
	}
	return count;
}
void sample_ribopos(vector<Ribopos>& sampled_ribopos, const vector<Ribopos>& ribopos, int count, const set<int>& included_columns)
{
	vector<int> all_hits;
	for (int i = 0; i < ribopos.size(); i++)
	{
		for (int j = 0; j < ribopos[i].reads[0]; j++)
		{
			all_hits.push_back(i);
		}
	}
	random_shuffle(all_hits.begin(), all_hits.end());
	if (count < all_hits.size())
	{
		all_hits = vector<int>(all_hits.begin(), all_hits.begin() + count);
	}
	sampled_ribopos = ribopos;
	for (int i = 0; i < ribopos.size(); i++)
	{
		sampled_ribopos[i].reads[0] = 0;
	}

	for (int i = 0; i < all_hits.size(); i++)
	{
		sampled_ribopos[all_hits[i]].reads[0]++;
	}
}



void identify_translated_orfs(string mod)
{

	vector<Experiment> experiments;
	//read_riboseq_studies(experiments, "riboseq_experiments_dataset.txt");
	//read_riboseq_studies(experiments, "/home/acwach/Riboseq/riboseq_experiments_dataset.txt");
	read_riboseq_studies(experiments, "input_files/riboseq_experiments_dataset.txt");
	//1.read ORFs
	vector<Orf> orfs;
	read_orfs(orfs, "orfs_comp", 1);

	cout << "\norfs read: " << orfs.size();
	//read ribosome profiling data
	vector<Ribopos> ribopos_f;
	vector<Ribopos> ribopos_r;

	vector<string> studies;
	read_riboseq_reads(ribopos_f, "input_files/mapped_ribseqs_combo_all_f", 0, -1, -1, studies);//data file on mapped riboseq reads on + strand 
	cout << "\npositions read: " << ribopos_f.size() << " " << studies.size();

	read_riboseq_reads(ribopos_r, "input_files/mapped_ribseqs_combo_all_r", 0, -1, -1, studies);//data file on mapped riboseq reads on - strand
	cout << "\nreverse positions read: " << ribopos_r.size() << " " << studies.size();

	if (mod=="multi_studies")//output riboseq read data mapped to ORFs grouped by riboseq study 
	{
		cout << "\nriboseq analysis: multi_studies";
		map<string, string> srp_map;
		set<string> studies_set;
		for (int i = 0; i < experiments.size(); i++)
		{
			studies_set.insert(experiments[i].srp);
			srp_map[experiments[i].srr] = experiments[i].srp;
		}
		for (auto it = studies_set.begin(); it != studies_set.end(); it++)
		{
			set<int> included_columns;
			for (int k = 0; k < studies.size(); k++)
			{
				if (srp_map.at(studies[k]) == *it)
				{
					included_columns.insert(k);
				}
			}
			map_riboseq_to_orfs(orfs, ribopos_f, 0, 1, included_columns);
			map_riboseq_to_orfs(orfs, ribopos_r, 1, 1, included_columns);
			cout << "\nmap riboseq to orfs";
			scramble_riboseq_reads(orfs);
			vector<riboseq_scorer> scorer;
			calc_riboseq_score(scorer, orfs);
			print_riboseqs_by_orf(scorer, orfs, "_" + (*it));
		}
	}
	else if (mod=="study_accumulation")//output riboseq read data mapped to ORFs grouped by number of studies considered, with rising number of studies. Used to assess how identification improves with more study effort
	{
		cout << "\nriboseq analysis: study_accumulation";
		map<string, string> srp_map;
		set<string> studies_set;
		vector<string> studies_vector;
		for (int i = 0; i < experiments.size(); i++)
		{
			studies_set.insert(experiments[i].srp);
			srp_map[experiments[i].srr] = experiments[i].srp;
		}
		for (auto it = studies_set.begin(); it != studies_set.end(); it++)
		{
			studies_vector.push_back(*it);
		}
		for(int a=0;a<50;a++)
		{
			random_shuffle(studies_vector.begin(),studies_vector.end());
			cout<<"\ntotal study count: "<<studies_vector.size();
			set<string> included_studies;
			for(int i=0;i<studies_vector.size();i++)
			{
				included_studies.insert(studies_vector[i]);
				set<int> included_columns;
				for (int k = 0; k < studies.size(); k++)
				{
					if (included_studies.count(srp_map.at(studies[k])))
					{
						included_columns.insert(k);
					}
				}
				map_riboseq_to_orfs(orfs, ribopos_f, 0, 1, included_columns);
				map_riboseq_to_orfs(orfs, ribopos_r, 1, 1, included_columns);
				cout << "\nmap riboseq to orfs";
				scramble_riboseq_reads(orfs);
				vector<riboseq_scorer> scorer;
				calc_riboseq_score(scorer, orfs);
				print_riboseqs_by_orf(scorer, orfs, "_accumulated_studies_" + to_string(a)+"_"+to_string(i)  );
			}
		}
	}
	else if (mod=="ypd_endpoint")//distinguish studies in YPD media vs other media
	{
		cout << "\nriboseq analysis: ypd_endpoint";
		map<string, int> ypd_map;
		for (int i = 0; i < experiments.size(); i++)
		{
			ypd_map[experiments[i].srr] = experiments[i].ypd;
		}
		vector<set<int>> included_columns(2);
		for (int s = 0; s < studies.size(); s++)
		{
			if (ypd_map[studies[s]] > -1)
			{
				if (ypd_map[studies[s]] == 1)
				{
					included_columns[1].insert(s);
				}
				else
				{
					included_columns[0].insert(s);
				}
			}
		}
		for (int i = 0; i < 2; i++)
		{
			vector<Orf> use_orfs = orfs;
			map_riboseq_to_orfs(use_orfs, ribopos_f, 0, 0, included_columns[i]);
			map_riboseq_to_orfs(use_orfs, ribopos_r, 1, 0, included_columns[i]);
			cout << "\nmap riboseq to orfs";
			scramble_riboseq_reads(use_orfs);
			vector<riboseq_scorer> scorer;
			calc_riboseq_score(scorer, use_orfs);
			print_riboseqs_by_orf(scorer, use_orfs, "ypd_endpoint" + to_string(i));
		}
	}
	else if (mod=="ypd_no_chx")//distinguish studies in YPD media with no CHX cell treatment
	{
		cout << "\nriboseq analysis: ypd_no_chx";
		map<string, int> ypd_map;
		map<string, int> chx_map;
		for (int i = 0; i < experiments.size(); i++)
		{
			ypd_map[experiments[i].srr] = experiments[i].ypd;
			chx_map[experiments[i].srr] = experiments[i].chx;
		}
		vector<set<int>> included_columns(2);
		for (int s = 0; s < studies.size(); s++)
		{
			if (ypd_map[studies[s]] > -1 && chx_map[studies[s]] > -1)
			{
				if (ypd_map[studies[s]] == 1 && chx_map[studies[s]] == 0)
				{
					included_columns[1].insert(s);
				}
				else
				{
					included_columns[0].insert(s);
				}
			}
		}
		for (int i = 0; i < 2; i++)
		{
			vector<Orf> use_orfs = orfs;
			map_riboseq_to_orfs(use_orfs, ribopos_f, 0, 0, included_columns[i]);
			map_riboseq_to_orfs(use_orfs, ribopos_r, 1, 0, included_columns[i]);
			cout << "\nmap riboseq to orfs";
			scramble_riboseq_reads(use_orfs);
			vector<riboseq_scorer> scorer;
			calc_riboseq_score(scorer, use_orfs);
			print_riboseqs_by_orf(scorer, use_orfs, "ypd_nochx_endpoint" + to_string(i));
		}
	}

	/*else if (chx_studies)
	{
		map<string, int> chx_map;
		for (int i = 0; i < experiments.size(); i++)
		{
			chx_map[experiments[i].srr] = experiments[i].chx;
		}
		vector<set<int>> included_columns(2);
		for (int s = 0; s < studies.size(); s++)
		{
			cout << "\n" << s << " " << studies[s] << " " << chx_map[studies[s]] << " " << studies.size();
			if (chx_map[studies[s]] > -1)
			{
				included_columns[chx_map[studies[s]]].insert(s);
			}
		}
		for (int i = 0; i < 2; i++)
		{
			vector<Orf> use_orfs = orfs;
			map_riboseq_to_orfs(use_orfs, ribopos_f, 0, 0, included_columns[i]);
			map_riboseq_to_orfs(use_orfs, ribopos_r, 1, 0, included_columns[i]);
			cout << "\nmap riboseq to orfs";
			scramble_riboseq_reads(use_orfs);
			vector<riboseq_scorer> scorer;
			calc_riboseq_score(scorer, use_orfs);
			print_riboseqs_by_orf(scorer, use_orfs, "chx" + to_string(i));
		}
	}*/
	else if (mod=="chx_studies_sampled")//distinguish studies with and without CHX treatment, sampling equal read depth for each category 
	{
		cout << "\nriboseq analysis: chx_studies_sampled";
		map<string, int> chx_map;
		for (int i = 0; i < experiments.size(); i++)
		{
			chx_map[experiments[i].srr] = experiments[i].chx;
		}
		vector<set<int>> included_columns(2);
		for (int s = 0; s < studies.size(); s++)
		{
			if (chx_map[studies[s]] > -1)
			{
				included_columns[chx_map[studies[s]]].insert(s);
			}
		}
		//cout << "\nnum exp:" << included_columns[0].size() << " " << included_columns[1].size();
		//getchar();
		vector<vector<Ribopos>> combined_ribopos_f(2);
		vector<vector<Ribopos>> combined_ribopos_r(2);
		vector<int> read_counts_f(2);
		vector<int> read_counts_r(2);
		set<int> column0;
		column0.insert(0);
		for (int i = 0; i < 2; i++)
		{
			assemble_ribopos(combined_ribopos_f[i], ribopos_f, included_columns[i]);
			assemble_ribopos(combined_ribopos_r[i], ribopos_r, included_columns[i]);
			read_counts_f[i] = count_reads(combined_ribopos_f[i]);
			read_counts_r[i] = count_reads(combined_ribopos_r[i]);
		}
		int sample_count_f = read_counts_f[0];
		int sample_count_r = read_counts_r[0];
		if (read_counts_f[1] < sample_count_f)
		{
			sample_count_f = read_counts_f[1];
		}
		if (read_counts_r[1] < sample_count_r)
		{
			sample_count_r = read_counts_r[1];
		}
		cout << "\nsample counts: " << sample_count_f << " " << sample_count_r << " " << read_counts_f[0] << "" << read_counts_f[1] << " " << included_columns[0].size() << " " << included_columns[1].size();
		for (int i = 0; i < 2; i++)
		{
			vector<Ribopos> sample_ribopos_f;
			vector<Ribopos> sample_ribopos_r;
			sample_ribopos(sample_ribopos_f, combined_ribopos_f[i], sample_count_f, included_columns[i]);
			sample_ribopos(sample_ribopos_r, combined_ribopos_r[i], sample_count_r, included_columns[i]);
			vector<Orf> use_orfs = orfs;
			map_riboseq_to_orfs(use_orfs, sample_ribopos_f, 0, 0, column0);// included_columns[i]);
			map_riboseq_to_orfs(use_orfs, sample_ribopos_r, 1, 0, column0);// included_columns[i]);
			cout << "\nmap riboseq to orfs";
			scramble_riboseq_reads(use_orfs);
			vector<riboseq_scorer> scorer;
			calc_riboseq_score(scorer, use_orfs);
			print_riboseqs_by_orf(scorer, use_orfs, "_chx_sampled" + to_string(i));
		}
	}
	else if (mod=="ypd_studies") //distinguish studies by media, sampling equal read depth for each category and considering a range of read depths
	{
		cout << "\nriboseq analysis: ypd_studies";
		vector<set<int>> included_columns(5);
		map<string, int> ypd_map;
		for (int i = 0; i < experiments.size(); i++)
		{
			ypd_map[experiments[i].srr] = experiments[i].ypd;
		}
		for (int s = 0; s < studies.size(); s++)
		{
			if (ypd_map[studies[s]] > -1)
			{
				included_columns[ypd_map[studies[s]]].insert(s);
				included_columns[2].insert(s);
			}
			included_columns[3].insert(s);
		}

		vector<vector<Ribopos>> combined_ribopos_f(5);
		vector<vector<Ribopos>> combined_ribopos_r(5);
		set<int> column0;
		column0.insert(0);
		for (int i = 0; i < 5; i++)
		{//random_ypd_by_experiment
			assemble_ribopos(combined_ribopos_f[i], ribopos_f, included_columns[i]);
			assemble_ribopos(combined_ribopos_r[i], ribopos_r, included_columns[i]);
			for (int s = 0; s < 20; s++)
			{
				cout << "\nis: " << i << " " << s << " " << ribopos_f.size();
				vector<Ribopos> sample_ribopos_f;
				vector<Ribopos> sample_ribopos_r;
				sample_ribopos(sample_ribopos_f, combined_ribopos_f[i], s * 10000000, included_columns[i]);
				sample_ribopos(sample_ribopos_r, combined_ribopos_r[i], s * 10000000, included_columns[i]);
				vector<Orf> use_orfs = orfs;
				map_riboseq_to_orfs(use_orfs, sample_ribopos_f, 0, 0, column0);
				map_riboseq_to_orfs(use_orfs, sample_ribopos_r, 1, 0, column0);
				cout << "\nmap riboseq to orfs";
				scramble_riboseq_reads(use_orfs);
				vector<riboseq_scorer> scorer;
				calc_riboseq_score(scorer, use_orfs);
				print_riboseqs_by_orf(scorer, use_orfs, "_ypd" + to_string(i) + "_" + to_string(s));
			}
		}
	}
	else //pool all studies togetether
	{
		set<int> included_columns;
		map_riboseq_to_orfs(orfs, ribopos_f, 0, -1, included_columns);
		map_riboseq_to_orfs(orfs, ribopos_r, 1, -1, included_columns);
		cout << "\nmap riboseq to orfs";
		scramble_riboseq_reads(orfs);
		vector<riboseq_scorer> scorer;

		calc_riboseq_score(scorer, orfs);
		print_riboseqs_by_orf(scorer, orfs, "");
	}
}

bool is_bad_syntenic_alignment(const Orf& orf, int species_id)
{
	if (orf.align_evalues[species_id] > .01)
	{
		return true;
	}
	int delcount = 0;
	for (int i = 0; i < orf.align_seq[species_id].size(); i++)
	{
		delcount += (orf.align_seq[species_id][i] == 4);
	}
	double del_freq = (double)delcount / (double)(1 + orf.end_pos - orf.start_pos);
	if (del_freq >= .4)
	{
		return true;
	}
	return false;
}

void get_blast_blocks(vector<MultiBlock>& blocks, const vector<Orf>& orfs, const vector<Blast_info>& blasts, int species_id, const vector<Contig>& contigs, const vector<Contig>& alt_contigs)
{
	map<string, int> contig_map;
	for (int i = 0; i < alt_contigs.size(); i++)
	{
		contig_map[alt_contigs[i].id] = i;
	}
	map<string, int> orf_map;
	for (int i = 0; i < orfs.size(); i++)
	{
		orf_map[orfs[i].orf_id] = i;
	}
	for (int i = 0; i < blasts.size(); i++)
	{
		if (!orf_map.count(blasts[i].str_id0))
		{
			cout << "\nbad map: " << blasts[i].str_id0;
			getchar();
		}
		int orf0 = orf_map.at(blasts[i].str_id0);
		if (is_bad_syntenic_alignment(orfs[orf0], species_id)/*trueorfs[i].good_align[species_id] == 0orfs[i].align_evalues[species_id]>1 && top_hit[i] != -1*/)
		{
			blocks.push_back(MultiBlock(2));
			blocks.back().specific_orf = orfs[orf0].orf_id;
			blocks.back().blast_id=i;
			blocks.back().blast_evalue=blasts[i].evalue;
			blocks.back().orf_id = orf0;
			blocks.back().start_gene = "X";
			blocks.back().end_gene = "X";
			blocks.back().start_pos[0] = orfs[orf0].start_pos - 1000;
			if (blocks.back().start_pos[0] < 0)
			{
				blocks.back().start_pos[0] = 0;
			}
			blocks.back().end_pos[0] = orfs[orf0].end_pos + 1000;
			if (blocks.back().end_pos[0] >= contigs[orfs[orf0].contig].seq.size())
			{
				blocks.back().end_pos[0] = contigs[orfs[orf0].contig].seq.size() - 1;
			}
			blocks.back().start_gene_start[0] = blocks.back().start_pos[0];
			blocks.back().start_gene_end[0] = orfs[orf0].start_pos - 1;
			blocks.back().end_gene_start[0] = orfs[orf0].end_pos + 1;
			blocks.back().end_gene_end[0] = blocks.back().end_pos[0];
			blocks.back().contig[0] = orfs[orf0].contig;
			blocks.back().start_pos[1] = blasts.at(i).align1_start - blasts.at(i).align0_start - 1000;//
			if (blocks.back().start_pos[1] < 0)
			{
				blocks.back().start_pos[1] = 0;
			}
			blocks.back().end_pos[1] = blasts.at(i).align1_end + (orfs[orf0].end_pos - orfs[orf0].start_pos - blasts.at(i).align0_end) + 1000;
			if (!contig_map.count(blasts.at(i).str_id1))
			{
				cout << "\nbad contig: " << blasts.at(i).str_id1;
				getchar();
			}
			blocks.back().contig[1] = contig_map.at(blasts.at(i).str_id1);
			if (blocks.back().end_pos[1] >= alt_contigs[blocks.back().contig[1]].seq.size())
			{
				blocks.back().end_pos[1] = alt_contigs[blocks.back().contig[1]].seq.size() - 1;
			}
			blocks.back().start_gene_start[1] = blocks.back().start_pos[1];
			blocks.back().start_gene_end[1] = blocks.back().start_pos[1] + 1000;
			blocks.back().end_gene_end[1] = blocks.back().end_pos[1];
			blocks.back().end_gene_start[1] = blocks.back().end_gene_end[1] - 1000;
			if (!contig_map.count(blasts.at(i).str_id1))
			{
				cout << "\n" << blasts.at(i).str_id1;
				getchar();
			}
			if (blasts.at(i).strand != orfs[orf0].strand)
			{
				blocks.back().rc = true;
			}
		}
	}
}

void read_blast_results(vector<Blast_info>& blasts, string filename, bool include_bitscore, bool include_lengths)
{
	ifstream file(filename);
	string line;
	while (getline(file, line))
	{
		vector<string> column_data;
		split(line, '\t', column_data);
		blasts.push_back(Blast_info());
		blasts.back().str_id0 = column_data[0];
		blasts.back().str_id1 = column_data[1];
		if (is_integer(column_data[0].substr(1)))
		{
			blasts.back().id0 = stoi(column_data[0].substr(1));
			blasts.back().id1 = stoi(column_data[1].substr(1));
		}
		blasts.back().identity = stod(column_data[2]);
		blasts.back().align_length = stoi(column_data[3]);
		blasts.back().mismatch_count = stoi(column_data[4]);
		blasts.back().gap_count = stoi(column_data[5]);
		blasts.back().align0_start = stoi(column_data[6]);
		blasts.back().align0_end = stoi(column_data[7]);
		blasts.back().align1_start = stoi(column_data[8]);
		blasts.back().align1_end = stoi(column_data[9]);
		if (blasts.back().align1_start > blasts.back().align1_end)
		{
			blasts.back().strand = 1;
			blasts.back().align1_start = stoi(column_data[9]);
			blasts.back().align1_end = stoi(column_data[8]);
		}
		else
		{
			blasts.back().strand = 0;
		}
		blasts.back().evalue = stod(column_data[10]);
		if (include_bitscore)
		{
			blasts.back().bitscore = stod(column_data[11]);
		}
		if (include_lengths)
		{
			blasts.back().length0 = stoi(column_data[11]);
			blasts.back().length1 = stoi(column_data[12]);
		}
	}
}

void filter_blasts(vector<Blast_info>& blasts, double threshold)
{
	vector<Blast_info> filtered_blasts;
	for (int i = 0; i < blasts.size(); i++)
	{
		if (blasts[i].evalue < threshold)
		{
			filtered_blasts.push_back(blasts[i]);
		}
		//.0000001
	}
	blasts = filtered_blasts;
}

void print_blocks(const vector<MultiBlock>& blocks, string filename)
{
	ofstream file(filename);
	file << "id start_gene end_gene specific_orf contig_focal start_focal end_focal rc blast_id blast_evalue start_gene_start_focal start_gene_end_focal end_gene_start_focal end_gene_end_focal";
	for (int i = 1; i < blocks[0].start_pos.size(); i++)
	{
		file << " contig" << i;
	}
	for (int i = 1; i < blocks[0].start_pos.size(); i++)
	{
		file << " start" << i;
	}
	for (int i = 1; i < blocks[0].start_pos.size(); i++)
	{
		file << " end" << i;
	}
	for (int i = 1; i < blocks[0].start_pos.size(); i++)
	{
		file << " start_gene_start" << i;
	}
	for (int i = 1; i < blocks[0].start_pos.size(); i++)
	{
		file << " start_gene_end" << i;
	}
	for (int i = 1; i < blocks[0].start_pos.size(); i++)
	{
		file << " end_gene_start" << i;
	}
	for (int i = 1; i < blocks[0].start_pos.size(); i++)
	{
		file << " end_gene_end" << i;
	}
	file << " seq_focal";
	for (int i = 0; i < blocks[0].start_pos.size(); i++)
	{
		file << " seq" << i;
	}
	for (int i = 0; i < blocks.size(); i++)
	{
		file << "\n" << i << " " << blocks[i].start_gene << " " << blocks[i].end_gene << " " << blocks[i].specific_orf << " " << blocks[i].contig[0] << " " << blocks[i].start_pos[0] << " " << blocks[i].end_pos[0] << " " << blocks[i].rc<<" "<<blocks[i].blast_id<<" "<<blocks[i].blast_evalue<<" "<< blocks[i].start_gene_start[0] << " " << blocks[i].start_gene_end[0] << " " << blocks[i].end_gene_start[0] << " " << blocks[i].end_gene_end[0];
		for (int j = 1; j < blocks[i].contig.size(); j++)
		{
			file << " " << blocks[i].contig[j];
		}
		for (int j = 1; j < blocks[i].start_pos.size(); j++)
		{
			file << " " << blocks[i].start_pos[j];
		}
		for (int j = 1; j < blocks[i].end_pos.size(); j++)
		{
			file << " " << blocks[i].end_pos[j];
		}
		for (int j = 1; j < blocks[i].start_gene_start.size(); j++)
		{
			file << " " << blocks[i].start_gene_start[j];
		}
		for (int j = 1; j < blocks[i].start_gene_end.size(); j++)
		{
			file << " " << blocks[i].start_gene_end[j];
		}
		for (int j = 1; j < blocks[i].end_gene_start.size(); j++)
		{
			file << " " << blocks[i].end_gene_start[j];
		}
		for (int j = 1; j < blocks[i].end_gene_end.size(); j++)
		{
			file << " " << blocks[i].end_gene_end[j];
		}
		for (int j = 0; j < blocks[i].align_seq.size(); j++)
		{
			file << " ";
			for (int k = 0; k < blocks[i].align_seq[j].size(); k++)
			{
				file << blocks[i].align_seq[j][k];
			}
		}
	}
}

void get_block_sequence(vector<MultiBlock>& blocks, const vector<Contig>& focal_contigs, const vector<Contig>& alt_contigs)
{
	for (int i = 0; i < blocks.size(); i++)
	{
		if(blocks[i].specific_orf=="0_22373_22429_0")
		{
			cout<<"\nblockmatch: "<<blocks[i].start_pos[1]<<" "<<blocks[i].end_pos[1]<<" "<<blocks[i].contig[1]<<" "<<blocks[i].rc;
			getchar();
		}

		for (int k = blocks[i].start_pos[0]; k <= blocks[i].end_pos[0]; k++)
		{
			blocks[i].seq[0].push_back(focal_contigs[blocks[i].contig[0]].seq.at(k));
		}
		for (int k = blocks[i].start_pos[1]; k <= blocks[i].end_pos[1]; k++)
		{
			blocks[i].seq[1].push_back(alt_contigs[blocks[i].contig[1]].seq.at(k));
		}
		if (blocks[i].rc)
		{
			reverse_complement(blocks[i].seq[1]);
		}
	}
}

void print_mult_fasta(string filename, const vector<vector<int>>& seq, const vector<string>& labels, bool print_gaps)
{
	map<int, char> rev_nuc_map = { { 0,'A' },{ 1,'C' },{ 2,'G' },{ 3,'T' },{ 4,'-' },{ 5,'N' },{ 6,'N' } };
	ofstream file(filename);
	for (int i = 0; i < seq.size(); i++)
	{
		file << ">" << labels[i] << "\n";
		for (int j = 0; j < seq[i].size(); j++)
		{
			if (!rev_nuc_map.count(seq[i][j]))
			{
				cout << "\nwrong nuc: " << seq[i][j] << " " << i << " " << filename << " " << labels[i];
				getchar();
			}
			else if (print_gaps || seq[i][j] != 4)
			{
				file << rev_nuc_map.at(seq[i][j]);
			}
		}
		file << "\n";
	}
}

bool read_muscle(vector<vector<int>>& seqs, map<string, int> labels, string filename)
{
	ifstream file(filename);
	map<char, int> nucmap = { { 'A',0 },{ 'a',0 },{ 'C',1 },{ 'c',1 },{ 'G',2 },{ 'g',2 },{ 'T',3 },{ 't',3 },{ 'N',5 },{ 'n',5 },{ '-',4 } };
	string line;
	int id = -1;
	while (getline(file, line))
	{
		if (line.substr(0, 1) == ">")
		{
			if (line.substr(2) == "f")
			{
				id = -1;
			}
			else
			{
				id = labels.at(line.substr(1));
			}
		}
		else
		{
			for (int i = 0; i < line.size(); i++)
			{
				if (id == -1)
				{
				}
				else if (!nucmap.count(line[i]))
				{
					return false;
				}
				else
				{
					seqs[id].push_back(nucmap.at(line[i]));
				}
			}
		}
	}
	return true;
}

void align_block_sequence(vector<MultiBlock>& blocks, string suffix, const vector<string>& labels, int from_align, int to_align)
{
	map<string, int> labels_map;
	for (int i = 0; i < labels.size(); i++)
	{
		labels_map[labels[i]] = i;
	}
	for (int i = from_align; i < to_align && i < blocks.size(); i++)
	{
		if (i > from_align&& blocks[i].start_pos[0] == blocks[i - 1].start_pos[0] && blocks[i].end_pos[0] == blocks[i - 1].end_pos[0] && blocks[i].contig[0] == blocks[i - 1].contig[0]
			&& blocks[i].start_pos[1] == blocks[i - 1].start_pos[1] && blocks[i].end_pos[1] == blocks[i - 1].end_pos[1] && blocks[i].contig[1] == blocks[i - 1].contig[1])
		{
			blocks[i].align_seq = blocks[i - 1].align_seq;
			blocks[i].align_success = blocks[i - 1].align_success;
		}
		else
		{
			print_mult_fasta("w" + to_string(from_align) + labels[1] + ".fasta", blocks[i].seq, labels, true);
			string cmd = "./muscle3.8.31_i86linux64 -in w" + to_string(from_align) + labels[1] + ".fasta -fastaout wm" + to_string(from_align) + labels[1] + ".out -quiet";
			int s = system(cmd.c_str());
			cout << "\nmult fasta printed";
			blocks[i].align_success = read_muscle(blocks[i].align_seq, labels_map, "wm" + to_string(from_align) + labels[1] + ".out");
			cout << "\nalign block: " << i << " " << blocks[i].align_seq[0].size() << " " << blocks[i].align_success << "\n" << cmd;
		}
	}
}

void blast_align(int species_id)
{
	vector<string> spp = { "Scer","Spar","Smik","Sjur","Skud","Sarb","Suva","Seub" };
	vector<string> genomes =
	{
		"S288C_reference_sequence_R64-2-1_20150113.fsa" ,
		"Spar.ultrascaf",
		"Smik.ultrascaf",
		"GCA_900290405.1_SacJureiUoM1_genomic.fna",
		"Skud.ultrascaf",
		"GCF_000292725.1_SacArb1.0_genomic.fna",
		"Sbay.ultrascaf",
		"GCF_001298625.1_SEUB3.0_genomic.fna"
	};

	vector<Orf> orfs;
	read_orfs_long(orfs, "orfs_scer_synteny_matched_checked", spp.size());
	cout << "\norfs read: " << orfs.size()<<" "<<orfs[0].orf_id;
	print_orfs(orfs, "_test");

	vector<vector<Contig>> contigs(spp.size());
	for (int i = 0; i < spp.size(); i++)
	{
		read_fasta(contigs[i], genomes[i]);
	}
	
	///
	///
	//cout<<"\nseqX:\n"; 
	//for(int q=0;q<40;q++)
	//{
	//	cout<<contigs[1][15].seq[590050+q];
	//}
	//getchar();
	///
	///
	
	vector<vector<MultiBlock>> blocks(spp.size());
	vector<vector<Blast_info>> blasts(spp.size());
	int i = species_id;
	cout << "\nspp: " << i;
	read_blast_results(blasts[i], "blasts_genome_nuc_" + spp[i] + "_3.out", false, false);
	cout << "\nblasts read: " << blasts[i].size();
	filter_blasts(blasts[i], .01);
	get_blast_blocks(blocks[i], orfs, blasts[i], i, contigs[0], contigs[i]);
	cout << "\nblocks constructed: " << blocks[i].size() << " out of " << orfs.size() << " orfs";
	print_blocks(blocks[i], "blocks_blast_" + to_string(i));
	get_block_sequence(blocks[species_id], contigs[0], contigs[species_id]);
	align_block_sequence(blocks[species_id], to_string(species_id), { spp[0],spp[species_id] }, 0, blocks[species_id].size());
	print_blocks(blocks[species_id], "aligned_blocks_blast_" + to_string(species_id));
}



void get_blast_map(map<int, int>& blast_map, const vector<Blast_info>& blasts)
{
	map<int, vector<int>> query_to_blast;
	for (int i = 0; i < blasts.size(); i++)
	{
		query_to_blast[blasts[i].id0].push_back(i);
	}
	for (auto it = query_to_blast.begin(); it != query_to_blast.end(); it++)
	{
		double max_bitscore = 0;
		int highest_score = 0;
		for (int i = 0; i < it->second.size(); i++)
		{
			if (blasts[(it->second).at(i)].bitscore > max_bitscore)
			{
				max_bitscore = blasts[(it->second).at(i)].bitscore;
				highest_score = i;
			}
		}
		blast_map[it->first] = blasts[(it->second).at(highest_score)].id1;
	}
}

void get_homology(vector<vector<Orf>>& orfs, const vector<map<int, int>>& blast_map, int focal_species)
{
	for (int i = 0; i < orfs[focal_species].size(); i++)
	{
		int homo_count = 0;
		for (int j = 0; j < blast_map.size(); j++)
		{
			if (j != focal_species && blast_map[j].count(i))
			{
				homo_count++;
			}
		}
		if (homo_count == orfs.size() - 1)
		{
			orfs[focal_species][i].homolog = true;
			for (int j = 0; j < orfs.size(); j++)
			{
				if (j != focal_species)
				{
					orfs[j][blast_map[j].at(i)].homolog = true;
					orfs[j][blast_map[j].at(i)].is_gene = orfs[focal_species][i].is_gene;
				}
			}
		}
	}
}

void filter_orfs_by_homology(vector<Orf>& orfs, map<int, int>& id_map)
{
	vector<Orf> filtered_orfs;
	for (int i = 0; i < orfs.size(); i++)
	{
		if (orfs[i].homolog)
		{
			id_map[orfs[i].id] = filtered_orfs.size();
			filtered_orfs.push_back(orfs[i]);
		}
	}
	orfs = filtered_orfs;
}

void get_orf_blocks(vector<MultiBlock>& blocks, const vector<Orf>& focal_orfs, const vector<Orf>& alt_orfs, const vector<Orf>& all_orfs, const map<int, int>& blast_map, const map<int, int>& id_map, const vector<Contig>& focal_contigs, const vector<Contig>& alt_contigs)
{
	map<int, int> contig_start;
	map<int, int> contig_end;
	for (int i = 0; i < focal_orfs.size(); i++)
	{
		if (i == 0 || focal_orfs[i].contig != focal_orfs[i - 1].contig)
		{
			contig_start[focal_orfs[i].contig] = i;

		}
		if (i == focal_orfs.size() - 1 || focal_orfs[i].contig != focal_orfs[i + 1].contig)
		{
			contig_end[focal_orfs[i].contig] = i;
		}
	}
	for (int i = 0; i < all_orfs.size(); i++)
	{
		int min_dist = 60000;
		int up_id = -1;
		int down_id = -1;
		bool rc = false;
		if (!contig_start.count(all_orfs[i].contig) || !contig_end.count(all_orfs[i].contig))
		{
			continue;
		}
		for (int j = contig_start.at(all_orfs[i].contig); j < contig_end.at(all_orfs[i].contig); j++)
		{
			for (int k = j + 1; k < contig_end.at(all_orfs[i].contig); k++)
			{
				if (focal_orfs[j].end_pos<all_orfs[i].start_pos &&
					focal_orfs[k].start_pos>all_orfs[i].end_pos&&
					blast_map.count(focal_orfs[j].id) &&
					blast_map.count(focal_orfs[k].id) &&
					id_map.count(blast_map.at(focal_orfs[j].id)) &&
					id_map.count(blast_map.at(focal_orfs[k].id)) &&
					alt_orfs[id_map.at(blast_map.at(focal_orfs[j].id))].contig == alt_orfs[id_map.at(blast_map.at(focal_orfs[k].id))].contig)
				{
					if (alt_orfs[id_map.at(blast_map.at(focal_orfs[j].id))].end_pos < alt_orfs[id_map.at(blast_map.at(focal_orfs[k].id))].start_pos &&
						(focal_orfs[k].end_pos - focal_orfs[j].start_pos) + (alt_orfs[id_map.at(blast_map.at(focal_orfs[k].id))].end_pos - alt_orfs[id_map.at(blast_map.at(focal_orfs[j].id))].start_pos) < min_dist)
					{
						min_dist = (focal_orfs[k].end_pos - focal_orfs[j].start_pos) + (alt_orfs[id_map.at(blast_map.at(focal_orfs[k].id))].end_pos - alt_orfs[id_map.at(blast_map.at(focal_orfs[j].id))].start_pos);
						up_id = j;
						down_id = k;
						rc = false;
					}
					else if (alt_orfs[id_map.at(blast_map.at(focal_orfs[k].id))].end_pos < alt_orfs[id_map.at(blast_map.at(focal_orfs[j].id))].start_pos &&
						(focal_orfs[k].end_pos - focal_orfs[j].start_pos) + (alt_orfs[id_map.at(blast_map.at(focal_orfs[j].id))].end_pos - alt_orfs[id_map.at(blast_map.at(focal_orfs[k].id))].start_pos) < min_dist)
					{
						min_dist = (focal_orfs[k].end_pos - focal_orfs[j].start_pos) + (alt_orfs[id_map.at(blast_map.at(focal_orfs[j].id))].end_pos - alt_orfs[id_map.at(blast_map.at(focal_orfs[k].id))].start_pos);
						up_id = j;
						down_id = k;
						rc = true;
					}

				}
			}
		}


		if (up_id > -1)
		{
			int up_id1 = id_map.at(blast_map.at(focal_orfs[up_id].id));
			int down_id1 = id_map.at(blast_map.at(focal_orfs[down_id].id));
			blocks.push_back(MultiBlock(2));

			blocks.back().rc = rc;
			if (blocks.back().rc)
			{
				up_id1 = id_map.at(blast_map.at(focal_orfs[down_id].id));
				down_id1 = id_map.at(blast_map.at(focal_orfs[up_id].id));
			}
			blocks.back().specific_orf = all_orfs[i].orf_id;
			blocks.back().orf_id = i;
			blocks.back().start_gene = focal_orfs[up_id].is_gene;
			blocks.back().end_gene = focal_orfs[down_id].is_gene;
			blocks.back().start_pos[0] = focal_orfs[up_id].start_pos;
			blocks.back().end_pos[0] = focal_orfs[down_id].end_pos;
			blocks.back().contig[0] = focal_orfs[up_id].contig;
			blocks.back().start_gene_start[0] = focal_orfs[up_id].start_pos;
			blocks.back().start_gene_end[0] = focal_orfs[up_id].end_pos;
			blocks.back().end_gene_end[0] = focal_orfs[down_id].end_pos;
			blocks.back().end_gene_start[0] = focal_orfs[down_id].start_pos;
			blocks.back().start_pos[1] = alt_orfs[up_id1].start_pos;
			blocks.back().end_pos[1] = alt_orfs[down_id1].end_pos;
			blocks.back().contig[1] = alt_orfs[up_id1].contig;
			blocks.back().start_gene_start[1] = alt_orfs[up_id1].start_pos;
			blocks.back().start_gene_end[1] = alt_orfs[up_id1].end_pos;
			blocks.back().end_gene_start[1] = alt_orfs[down_id1].start_pos;
			blocks.back().end_gene_end[1] = alt_orfs[down_id1].end_pos;
		}
		else //try to anchor to end of contig
		{
			int min_up_dist = 80000;
			int min_down_dist = 80000;
			for (int j = contig_start.at(all_orfs[i].contig); j < contig_end.at(all_orfs[i].contig); j++)
			{
				if (focal_orfs[j].end_pos < all_orfs[i].start_pos &&
					blast_map.count(focal_orfs[j].id) &&
					id_map.count(blast_map.at(focal_orfs[j].id)) &&
					(focal_contigs[all_orfs[i].contig].seq.size() - 1 - focal_orfs[j].start_pos) + (alt_contigs[alt_orfs[id_map.at(blast_map.at(focal_orfs[j].id))].contig].seq.size() - 1 - alt_orfs[id_map.at(blast_map.at(focal_orfs[j].id))].start_pos) < min_up_dist
					)
				{
					min_up_dist = (focal_orfs[all_orfs[i].contig].seq.size() - 1 - focal_orfs[j].start_pos) + (alt_contigs[alt_orfs[id_map.at(blast_map.at(focal_orfs[j].id))].contig].seq.size() - 1 - alt_orfs[id_map.at(blast_map.at(focal_orfs[j].id))].start_pos);
					up_id = j;// orfs[0][j].id;
				}
				else if (focal_orfs[j].start_pos > all_orfs[i].end_pos&&
					blast_map.count(focal_orfs[j].id) &&
					id_map.count(blast_map.at(focal_orfs[j].id)) &&
					(focal_orfs[j].end_pos) + (alt_orfs[id_map.at(blast_map.at(focal_orfs[j].id))].end_pos) < min_down_dist
					)
				{
					min_down_dist = (focal_orfs[j].end_pos) + (alt_orfs[id_map.at(blast_map.at(focal_orfs[j].id))].end_pos);
					down_id = j;// orfs[0][j].id;
				}
			}
			if (up_id != -1)
			{
				int up_id1 = id_map.at(blast_map.at(focal_orfs[up_id].id));
				blocks.push_back(MultiBlock(2));
				blocks.back().specific_orf = all_orfs[i].orf_id;
				blocks.back().orf_id = i;
				blocks.back().start_gene = focal_orfs[up_id].is_gene;
				blocks.back().end_gene = "CONTIG_END";
				blocks.back().start_pos[0] = focal_orfs[up_id].start_pos;
				blocks.back().end_pos[0] = focal_contigs[focal_orfs[up_id].contig].seq.size() - 1;
				blocks.back().contig[0] = focal_orfs[up_id].contig;
				blocks.back().start_gene_start[0] = focal_orfs[up_id].start_pos;
				blocks.back().start_gene_end[0] = focal_orfs[up_id].end_pos;
				blocks.back().end_gene_end[0] = focal_contigs[focal_orfs[up_id].contig].seq.size() - 1;
				blocks.back().end_gene_start[0] = focal_contigs[focal_orfs[up_id].contig].seq.size() - 1;

				blocks.back().start_pos[1] = alt_orfs[up_id1].start_pos;
				blocks.back().end_pos[1] = alt_contigs[alt_orfs[up_id1].contig].seq.size() - 1;
				blocks.back().contig[1] = alt_orfs[up_id1].contig;
				blocks.back().start_gene_start[1] = alt_orfs[up_id1].start_pos;
				blocks.back().start_gene_end[1] = alt_orfs[up_id1].end_pos;
				blocks.back().end_gene_start[1] = alt_contigs[alt_orfs[up_id1].contig].seq.size() - 1;
				blocks.back().end_gene_end[1] = alt_contigs[alt_orfs[up_id1].contig].seq.size() - 1;
			}
			if (down_id != -1)
			{
				int down_id1 = id_map.at(blast_map.at(focal_orfs[down_id].id));
				blocks.push_back(MultiBlock(2));
				blocks.back().specific_orf = all_orfs[i].orf_id;
				blocks.back().orf_id = i;
				blocks.back().start_gene = "CONTIG_START";
				blocks.back().end_gene = focal_orfs[down_id].is_gene;
				blocks.back().start_pos[0] = 0;
				blocks.back().end_pos[0] = focal_orfs[down_id].end_pos;
				blocks.back().contig[0] = focal_orfs[down_id].contig;
				blocks.back().start_gene_start[0] = 0;
				blocks.back().start_gene_end[0] = 0;
				blocks.back().end_gene_end[0] = focal_orfs[down_id].end_pos;
				blocks.back().end_gene_start[0] = focal_orfs[down_id].start_pos;

				blocks.back().start_pos[1] = 0;
				blocks.back().end_pos[1] = alt_orfs[down_id1].end_pos;
				blocks.back().contig[1] = alt_orfs[down_id1].contig;
				blocks.back().start_gene_start[1] = 0;
				blocks.back().start_gene_end[1] = 0;
				blocks.back().end_gene_start[1] = alt_orfs[down_id1].start_pos;
				blocks.back().end_gene_end[1] = alt_orfs[down_id1].end_pos;
			}
		}
	}
}

void print_block_sequence_fasta(vector<MultiBlock>& blocks, string spp, string suffix)
{
	map<int, char> rev_nuc_map = { { 0,'A' },{ 1,'C' },{ 2,'G' },{ 3,'T' },{ 5,'N' },{ 6,'N' } };
	ofstream file("blocks_Scer_" + spp + suffix);
	for (int i = 0; i < blocks.size(); i++)
	{
		file << ">" << blocks[i].specific_orf << "," << blocks[i].start_gene << "," << blocks[i].end_gene << ",Scer\n";// << spp;
		for (int j = 0; j < blocks[i].seq[0].size(); j++)
		{
			file << rev_nuc_map.at(blocks[i].seq[0][j]);
		}
		file << "\n>" << blocks[i].specific_orf << "," << blocks[i].start_gene << "," << blocks[i].end_gene << "," << spp << "\n";
		for (int j = 0; j < blocks[i].seq[1].size(); j++)
		{
			file << rev_nuc_map.at(blocks[i].seq[1][j]);
		}
		file << "\n";
	}
}

void synteny_align(int species_id, bool FORCE_MAKE_BLAST_DBS, bool FORCE_REDO_BLASTS)
{
	int focal_species = 0;
	vector<string> spp = { "Scer","Spar","Smik","Sjur","Skud","Sarb","Suva","Seub" };
	vector<string> genomes = {
		"S288C_reference_sequence_R64-2-1_20150113.fsa" ,
		"Spar.ultrascaf",
		"Smik.ultrascaf",
		"GCA_900290405.1_SacJureiUoM1_genomic.fna",
		"Skud.ultrascaf",
		"GCF_000292725.1_SacArb1.0_genomic.fna",
		"Sbay.ultrascaf",
		"GCF_001298625.1_SEUB3.0_genomic.fna"
	};
	vector<vector<Contig>> contigs(spp.size());
	vector<vector<Orf>> orfs(spp.size());
	for (int i = 0; i < spp.size(); i++)
	{
		read_fasta(contigs[i], genomes[i]);
		get_orfs(orfs[i], contigs[i], spp.size(), 200, true);
		string filename = "orfs_aa_" + spp[i] + ".fas";
		//2.a. setup blast databases if don't already exist
		if (!file_exists(filename) || FORCE_MAKE_BLAST_DBS)
		{
			print_orfs_aa_fasta(filename, orfs[i]);
			string cmd = "makeblastdb -in " + filename + " -dbtype prot -out " + spp[i];
			int s = system(cmd.c_str());
		}
	}

	map<string, int> contig_map;
	for (int i = 0; i < contigs[0].size(); i++)
	{
		contig_map[contigs[0][i].id] = i;
	}

	string aa_filename = "orfs_seqs" + to_string(focal_species) + ".fas";
	print_orfs_aa_fasta(aa_filename, orfs[focal_species]);
	for (int i = 0; i < spp.size(); i++)
	{
		print_orfs(orfs[i], spp[i]);
	}
	vector<vector<Blast_info>> blasts(spp.size());
	vector<map<int, int>> blast_map(spp.size());
	for (int i = 0; i < spp.size(); i++)
	{
		if (i == focal_species)
		{
			continue;
		}
		string filename = "blasts_" + spp[i] + "_" + to_string(focal_species) + ".out";
		if (!file_exists(filename) || FORCE_REDO_BLASTS)
		{
			string cmd = "blastp -db " + spp[i] + " -query " + aa_filename + " -out " + filename + " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen\"";
			cout << "\nrunning command: " << cmd;
			int s = system(cmd.c_str());
		}

		read_blast_results(blasts[i], filename, true, false);
		filter_blasts(blasts[i], .0000001);
		get_blast_map(blast_map[i], blasts[i]);
		cout << "\nblasts read for species " << i << ":" << blasts[i].size();
	}

	vector<Orf> genes;
	read_annotated_genes(genes, "orf_coding.fasta", true);
	
	check_orf_overlap(orfs[0], genes);

	//4. for each ORF in focal species, build syntenic block around it by finding upstream and downstream anchors in common among relatives
	vector<map<int, int>> id_map(spp.size());
	get_homology(orfs, blast_map, focal_species);
	for (int i = 0; i < orfs.size(); i++)
	{
		filter_orfs_by_homology(orfs[i], id_map[i]);
	}
	cout << "\nhomologous orfs: " << orfs[focal_species].size();

	vector<MultiBlock> blocks;
	vector<Orf> scer_orfs;
	get_orfs(scer_orfs, contigs[focal_species], spp.size(), 7, false);
	remove_duplicate_orfs_by_contig(scer_orfs, contigs[focal_species]);
	cout << "\n_orfs remaining after remove duplicates: " << scer_orfs.size();

	get_orf_blocks(blocks, orfs[focal_species], orfs[species_id], scer_orfs, blast_map[species_id], id_map[species_id], contigs[focal_species], contigs[species_id]);
	cout << "\nblocks identified: " << blocks.size();
	print_blocks(blocks, "blocks" + to_string(species_id) + "_" + to_string(focal_species));
	get_block_sequence(blocks, contigs[focal_species], contigs[species_id]);
	print_block_sequence_fasta(blocks, spp[species_id], "_" + to_string(focal_species));
	align_block_sequence(blocks, to_string(species_id), { spp[focal_species],spp[species_id] }, 0, blocks.size());
	print_blocks(blocks, "aligned_blocks" + to_string(species_id) + "_" + to_string(focal_species));
}



void read_multiblock(vector<MultiBlock>& blocks, string filename, int species_count)
{
	ifstream file(filename);
	string line;
	getline(file, line);
	while (getline(file, line))
	{
		int x = 1;
		vector<string> column_data;
		split(line, ' ', column_data);
		blocks.push_back(MultiBlock(species_count));
		blocks.back().start_gene = column_data[x++];
		blocks.back().end_gene = column_data[x++];
		blocks.back().specific_orf = column_data[x++];
		blocks.back().contig[0] = stoi(column_data[x++]);
		blocks.back().start_pos[0] = stoi(column_data[x++]);
		blocks.back().end_pos[0] = stoi(column_data[x++]);
		blocks.back().rc = stoi(column_data[x++]);
		blocks.back().blast_id = stoi(column_data[x++]);
		blocks.back().blast_evalue = stod(column_data[x++]);
		blocks.back().start_gene_start[0] = stoi(column_data[x++]);
		blocks.back().start_gene_end[0] = stoi(column_data[x++]);
		blocks.back().end_gene_start[0] = stoi(column_data[x++]);
		blocks.back().end_gene_end[0] = stoi(column_data[x++]);
		for (int i = 1; i < blocks.back().contig.size(); i++)
		{
			blocks.back().contig[i] = stoi(column_data[x++]);
		}
		for (int i = 1; i < blocks.back().start_pos.size(); i++)
		{
			blocks.back().start_pos[i] = stoi(column_data[x++]);
		}
		for (int i = 1; i < blocks.back().end_pos.size(); i++)
		{
			blocks.back().end_pos[i] = stoi(column_data[x++]);
		}
		for (int i = 1; i < blocks.back().start_gene_start.size(); i++)
		{
			blocks.back().start_gene_start[i] = stoi(column_data[x++]);
		}
		for (int i = 1; i < blocks.back().start_gene_end.size(); i++)
		{
			blocks.back().start_gene_end[i] = stoi(column_data[x++]);
		}
		for (int i = 1; i < blocks.back().end_gene_start.size(); i++)
		{
			blocks.back().end_gene_start[i] = stoi(column_data[x++]);
		}
		for (int i = 1; i < blocks.back().end_gene_end.size(); i++)
		{
			blocks.back().end_gene_end[i] = stoi(column_data[x++]);
		}
		for (int i = 0; x < column_data.size() && i < column_data[x].size(); i++)
		{
			blocks.back().align_seq[0].push_back(column_data[x][i] - '0');
		}
		x++;
		for (int s = 1; s < blocks.back().align_seq.size(); s++)
		{
			for (int i = 0; x < column_data.size() && i < column_data[x].size(); i++)
			{
				blocks.back().align_seq[s].push_back(column_data[x][i] - '0');
			}
			x++;
		}
	}
}

void get_block_limits(MultiBlock& block)
{
	for (int i = 0; i < block.start_pos.size(); i++)
	{
		int index = -1;
		for (int j = 0; j < block.align_seq[i].size(); j++)
		{
			if (block.align_seq[i][j] != 4)
			{
				index++;
				if (index == block.start_gene_end.at(i) - block.start_pos[i])
				{
					block.in_start_intergene[i] = j;
				}
				if (index == block.end_gene_start[i] - block.start_pos[i])
				{
					block.in_end_intergene[i] = j;
				}
			}
		}
	}
}

void check_blocks_good(vector<MultiBlock>& blocks)
{
	for (int i = 0; i < blocks.size(); i++)
	{
		for (int j = 0; j < blocks[i].start_pos.size(); j++)
		{
			if (blocks[i].start_pos[j] > blocks[i].end_pos[j] && blocks[i].end_pos[0] - blocks[i].start_pos[0] < 50000)
			{
				blocks[i].is_good = false;
			}
			if (blocks[i].align_seq[0].size() != blocks[i].align_seq[1].size())
			{
				blocks[i].is_good = false;
			}
		}
	}
}

void revc_blocks(vector<MultiBlock>& blocks)
{
	for (int i = 0; i < blocks.size(); i++)
	{
		for (int s = 0; s < blocks[i].align_seq.size(); s++)
		{
			reverse_complement(blocks[i].align_seq[s]);
		}
		for (int s = 0; s < blocks[i].in_end_intergene.size(); s++)
		{
			int in_end = blocks[i].in_end_intergene[s];
			blocks[i].in_end_intergene[s] = blocks[i].align_seq[s].size() - blocks[i].in_start_intergene[s] - 1;
			blocks[i].in_start_intergene[s] = blocks[i].align_seq[s].size() - in_end - 1;
		}
	}
}

void read_sub_rates(vector<vector<double>>& sub_rates, string filename)
{
	ifstream file(filename);
	string line;
	getline(file, line);
	while (getline(file, line))
	{
		vector<string> rows;
		split(line, ' ', rows);
		sub_rates.push_back(vector<double>());
		for (int i = 1; i < rows.size(); i++)
		{
			sub_rates.back().push_back(stod(rows[i]));
		}
	}
}

void find_orf_in_block(const MultiBlock& block, Orf& orf, int strand, int species)
{
	if (orf.strand == strand || strand == -1)
	{
		int index = -1;
		int in_orf_position = -1;
		int alt_index = -1;
		int orf_start = orf.start_pos - block.start_pos[0];
		int orf_end = orf.end_pos - block.start_pos[0];
		if (orf.strand == 1 && strand != -1)
		{
			orf_start = block.end_pos[0] - orf.end_pos;
			orf_end = block.end_pos[0] - orf.start_pos;
		}
		int pos = 0;
		bool in_orf = false;
		for (int i = 0; i < block.align_seq[0].size(); i++)
		{
			if (block.align_seq.at(1).at(i) != 4)
			{
				alt_index++;
			}
			if (block.align_seq[0][i] != 4)
			{
				index++;
				if (index == orf_start)
				{
					orf.in_start = i;
					pos = 0;
					in_orf = true;
					orf.alt_contig[species] = block.contig[1];
					//if (strand == 0)//change to block strand
					if(block.rc==0)
					{
						orf.alt_start_pos[species] = alt_index + block.start_pos[1];
					}
					else
					{
						orf.alt_end_pos[species] = block.end_pos[1] - alt_index;
					}
				}
				if (in_orf)
				{
					in_orf_position++;
					orf.frame_seq.push_back(pos);
					orf.align_index.push_back(in_orf_position);
					pos++;
					if (pos == 3)
					{
						pos = 0;
					}
				}
				else
				{
					orf.frame_seq.push_back(5);
					orf.align_index.push_back(-1);
				}
				if (index == orf_end)
				{
					orf.in_end = i;
					in_orf = false;
					//if (strand == 0)
					if(block.rc==0)
					{
						orf.alt_end_pos[species] = alt_index + block.start_pos[1];
					}
					else
					{
						orf.alt_start_pos[species] = block.end_pos[1] - alt_index;
					}
					//return;
				}
			}

			else if (in_orf)
			{
				orf.frame_seq.push_back(4);
				orf.align_index.push_back(-1);
			}
			else
			{
				orf.frame_seq.push_back(5);
				orf.align_index.push_back(-1);
			}
		}
	}
}

void get_ref_pos(vector<int>& ref_pos, const vector<int>& align)
{
	ref_pos = vector<int>(align.size());
	for (int i = 1; i < align.size(); i++)
	{
		if (align[i] != 4)
		{
			ref_pos[i] = ref_pos[i - 1] + 1;
		}
		else
		{
			ref_pos[i] = ref_pos[i - 1];
		}
	}
}

bool codon_match(const vector<int>& test, const vector<vector<int>>& match_set)
{
	for (int i = 0; i < match_set.size(); i++)
	{
		bool all_good = true;
		for (int j = 0; j < match_set[i].size(); j++)
		{
			if (test[j] != match_set[i][j])
			{
				all_good = false;
				break;
			}
		}
		if (all_good)
		{
			return true;
		}
	}
	return false;
}

void aquire_orfs(vector<Orf>& orfs, const MultiBlock& block, int strand, const vector<int>& ref_pos)
{
	vector<vector<int>> start_codons = { {0,3,2} };
	vector<vector<int>> stop_codons = { {3,0,0},{3,0,2},{3,2,0} };
	int nuc0 = -1;
	int nuc1 = -1;
	int nuc2 = -1;
	int nuc0_pos = -1;
	int nuc1_pos = -1;
	int nuc2_pos = -1;
	vector<int> real_pos;
	real_pos.push_back(0);
	for (int j = 1; j < block.align_seq[1].size(); j++)
	{
		if (block.align_seq[1][j] != 4)
		{
			real_pos.push_back(real_pos.back() + 1);
		}
		else
		{
			real_pos.push_back(real_pos.back());
		}
	}
	for (int j = 0; j < block.align_seq[1].size(); j++)
	{
		if (block.align_seq[1][j] != 4)
		{
			nuc0 = nuc1;
			nuc1 = nuc2;
			nuc2 = block.align_seq[1][j];
			nuc0_pos = nuc1_pos;
			nuc1_pos = nuc2_pos;
			nuc2_pos = j;
			vector<int> codon = { nuc0,nuc1,nuc2 };
			if (codon_match(codon, start_codons))
			{
				Orf orf(block.start_pos.size());
				int nuc_pos0 = -1;
				int nuc_pos1 = -1;
				int nuc_pos2 = -1;
				int pos = 0;
				int cur_offset = 0;
				for (int k = nuc0_pos; k < block.align_seq[1].size(); k++)
				{
					if (block.align_seq[1][k] != 4)
					{
						orf.seq.push_back(block.align_seq[1][k]);
						nuc_pos0 = nuc_pos1;
						nuc_pos1 = nuc_pos2;
						nuc_pos2 = block.align_seq[1][k];
						vector<int> orf_codon = { nuc_pos0,nuc_pos1,nuc_pos2 };
						if (pos == 2 && codon_match(orf_codon, stop_codons))
						{
							orf.contig = block.contig[1];
							orf.strand = strand;
							orf.start_pos = real_pos[nuc0_pos] + block.start_pos[1];
							orf.end_pos = real_pos[k] + block.start_pos[1];
							orf.ref_start = ref_pos[nuc0_pos] + block.start_pos[0];
							orf.ref_end = ref_pos[k] + block.start_pos[0];
							orf.ref_contig = block.contig[0];
							orf.in_start = nuc0_pos;
							orf.in_end = k;
							if (strand == 1)
							{
								orf.end_pos = block.end_pos[1] - real_pos[nuc0_pos];
								orf.start_pos = block.end_pos[1] - real_pos[k];
								orf.ref_end = block.end_pos[0] - ref_pos[nuc0_pos];
								orf.ref_start = block.end_pos[0] - ref_pos[k];
							}
							orfs.push_back(orf);
							break;
						}
						pos++;
						if (pos == 3)
						{
							pos = 0;
						}
					}
					else
					{
					}
				}
			}
		}
	}
}

void find_alt_orf_in_block(const MultiBlock& block, Orf& orf, int strand)
{
	if (orf.strand == strand)
	{
		int pos = 0;
		bool in_orf = false;
		for (int i = 0; i < block.align_seq[1].size(); i++)
		{
			if (i == orf.in_start)
			{
				in_orf = true;
				pos = 0;
			}
			if (block.align_seq[1][i] != 4)
			{
				if (in_orf)
				{
					orf.frame_seq.push_back(pos);
					pos++;
					if (pos == 3)
					{
						pos = 0;
					}
				}
				else
				{
					orf.frame_seq.push_back(5);
				}
				if (i == orf.in_end)
				{
					in_orf = false;
				}
			}
			else if (in_orf)
			{
				orf.frame_seq.push_back(4);
			}
			else
			{
				orf.frame_seq.push_back(5);
			}
		}
	}
}

int count_N(const vector<int>& seq)
{
	int count = 0;
	for (int i = 0; i < seq.size(); i++)
	{
		if (seq[i] == 5)
		{
			count++;
		}
	}
	return count;
}

bool single_compare_orfs(Orf& focal_orf, MultiBlock& block, int species, const vector<vector<double>>& sub_rates)
{
	bool changed = false;
	map<int, char> aa_map = { { 0,'K' },{ 1,'N' },{ 2,'K' },{ 3,'N' },{ 4,'T' },{ 5,'T' },{ 6,'T' },{ 7,'T' },{ 8,'R' },{ 9,'S' },{ 10,'R' },{ 11,'S' },{ 12,'I' },{ 13,'I' },{ 14,'M' },{ 15,'I' },{ 16,'Q' },{ 17,'H' },{ 18,'Q' },{ 19,'H' },{ 20,'P' },{ 21,'P' },{ 22,'P' },{ 23,'P' },{ 24,'R' },{ 25,'R' },{ 26,'R' },{ 27,'R' },{ 28,'L' },{ 29,'L' },{ 30,'L' },{ 31,'L' },{ 32,'E' },{ 33,'D' },{ 34,'E' },{ 35,'D' },{ 36,'A' },{ 37,'A' },{ 38,'A' },{ 39,'A' },{ 40,'G' },{ 41,'G' },{ 42,'G' },{ 43,'G' },{ 44,'V' },{ 45,'V' },{ 46,'V' },{ 47,'V' },{ 48,'X' },{ 49,'Y' },{ 50,'X' },{ 51,'Y' },{ 52,'S' },{ 53,'S' },{ 54,'S' },{ 55,'S' },{ 56,'X' },{ 57,'C' },{ 58,'W' },{ 59,'C' },{ 60,'L' },{ 61,'F' },{ 62,'L' },{ 63,'F' } };
	if (block.is_good)
	{
		vector<Orf> orfs;
		vector<int> ref_pos;
		get_ref_pos(ref_pos, block.align_seq[0]);
		aquire_orfs(orfs, block, focal_orf.strand, ref_pos);
		//cout << "\norfs aq: " << orfs.size()<<" "<<block.align_seq[0].size();
		///getchar();
		//filter_orfs_by_same_stop(orfs);
		filter_orfs_by_length(orfs, 30);

		vector<int> focal_align_seq;
		vector<int> align_seq;
		for (int p = 0; p < focal_orf.frame_seq.size(); p++)
		{
			if (focal_orf.frame_seq[p] < 5)
			{
				focal_align_seq.push_back(block.align_seq.at(0).at(p));
				align_seq.push_back(block.align_seq.at(1).at(p));
			}
		}

		int gap_count = 0;
		int n_count = 0;
		int diff_count = 0;
		for (int p = 0; p < focal_align_seq.size(); p++)
		{
			if (focal_align_seq[p] != align_seq[p] && (focal_align_seq[p] == 4 || align_seq[p] == 4))
			{
				gap_count++;
			}
			if (focal_align_seq[p] != align_seq[p] && (focal_align_seq[p] == 5 || align_seq[p] == 5))
			{
				n_count++;
			}
			else if (focal_align_seq[p] != align_seq[p])
			{
				diff_count++;
			}
		}
		if (gap_count > .4 * (double)focal_align_seq.size() || n_count > 0 || focal_align_seq.size() == 0 || diff_count > .4 * focal_align_seq.size())
		{
			focal_orf.good_align[species] = 0;
		}
		else
		{
			focal_orf.good_align[species] = 1;
		}
		focal_orf.focal_align_seq[species] = focal_align_seq;
		focal_orf.align_seq[species] = align_seq;
		focal_orf.align_identity[species] = 1.0 - (double)diff_count / (double)focal_orf.focal_align_seq[species].size();

		for (int k = 0; k < orfs.size(); k++)
		{
			int max_overlap = get_max_overlap(focal_orf.in_start, focal_orf.in_end, orfs[k].in_start, orfs[k].in_end);
			if (max_overlap > focal_orf.overlap[species])
			{
				//check_overlap++;
				find_alt_orf_in_block(block, orfs[k], focal_orf.strand);
				vector<int> match_frame;
				vector<int> diffs(3);
				int blast_length = 0;
				//double evalue = -1;
				//get_blast_comparison(orfs[focal][j], orfs[s][k],evalue,blast_length);
				int overlap = 0;
				int syn_diffs = 0;
				int nonsyn_diffs = 0;
				double syn_exp = 0;
				double nonsyn_exp = 0;
				int total_muts = 0;
				double total_exp = 0;
				for (int p = 0; p < focal_orf.frame_seq.size(); p++)
				{
					if (focal_orf.frame_seq[p] < 4 && focal_orf.frame_seq[p] == orfs[k].frame_seq[p])
					{
						overlap++;
						if (block.align_seq[0][p] != block.align_seq[1][p])
						{
							diffs.at(focal_orf.frame_seq[p])++;
						}
					}
					if (focal_orf.frame_seq[p] < 5)
					{
						if (focal_orf.frame_seq[p] == 4 || orfs[k].frame_seq[p] == 4)
						{
							match_frame.push_back(4);
						}
						else if (focal_orf.frame_seq[p] == orfs[k].frame_seq[p])
						{
							match_frame.push_back(1);
							int nuc0 = block.align_seq[0][p];
							int nuc1 = block.align_seq[1][p];
							for (int q = 0; q < 4; q++)
							{
								if (q == nuc0)
								{
									continue;
								}
								int old_codon_hash = -1;
								int new_codon_hash = -1;
								int align_pos = focal_orf.align_index.at(p);
								if (focal_orf.frame_seq[p] == 0)
								{
									old_codon_hash = focal_orf.seq.at(align_pos) * 4 * 4 + focal_orf.seq.at(align_pos + 1) * 4 + focal_orf.seq.at(align_pos + 2);
									new_codon_hash = q * 4 * 4 + focal_orf.seq.at(align_pos + 1) * 4 + focal_orf.seq.at(align_pos + 2);
								}
								if (focal_orf.frame_seq[p] == 1)
								{
									old_codon_hash = focal_orf.seq.at(align_pos - 1) * 4 * 4 + focal_orf.seq.at(align_pos) * 4 + focal_orf.seq.at(align_pos + 1);
									new_codon_hash = focal_orf.seq.at(align_pos - 1) * 4 * 4 + q * 4 + focal_orf.seq.at(align_pos + 1);
								}
								if (focal_orf.frame_seq[p] == 2)
								{
									old_codon_hash = focal_orf.seq.at(align_pos - 2) * 4 * 4 + focal_orf.seq.at(align_pos - 1) * 4 + focal_orf.seq.at(align_pos);
									new_codon_hash = focal_orf.seq.at(align_pos - 2) * 4 * 4 + focal_orf.seq.at(align_pos - 1) * 4 + q;
								}
								if (aa_map.at(old_codon_hash) == aa_map.at(new_codon_hash))
								{
									if (nuc0 < 4 && nuc1 < 4)
									{
										syn_exp += sub_rates[nuc0][q];
										total_exp += sub_rates[nuc0][q];
										if (nuc1 == q)
										{
											syn_diffs++;
											total_muts++;
										}
									}
								}
								else
								{
									if (nuc0 < 4 && nuc1 < 4)
									{
										nonsyn_exp += sub_rates[nuc0][q];
										total_exp += sub_rates[nuc0][q];
										if (nuc1 == q)
										{
											nonsyn_diffs++;
											total_muts++;
										}
									}
								}
							}
						}
						else
						{
							match_frame.push_back(0);
						}
					}
				}

				if (overlap > focal_orf.overlap[species])
				{
					changed = true;
					focal_orf.overlap[species] = overlap;
					focal_orf.match_frame[species] = match_frame;
					focal_orf.diffs[species] = diffs;
					focal_orf.homo_length[species] = orfs[k].end_pos - orfs[k].start_pos;
					focal_orf.syn_diffs[species] = syn_diffs;
					focal_orf.nonsyn_diffs[species] = nonsyn_diffs;
					if (total_exp == 0)
					{
						focal_orf.syn_diffs_exp[species] = 0;
						focal_orf.nonsyn_diffs_exp[species] = 0;
					}
					else
					{
						focal_orf.syn_diffs_exp[species] = syn_exp * (double)total_muts / total_exp;
						focal_orf.nonsyn_diffs_exp[species] = nonsyn_exp * (double)total_muts / total_exp;
					}
					focal_orf.alt_seqs[species] = orfs[k].seq;
					if (count_N(orfs[k].seq) == 0)
					{
						focal_orf.alt_aaseq[species] = translate(orfs[k].seq);
					}
				}
			}
		}
	}
	return changed;
}

void read_orf_bounds(vector<Orf>& orfs, string filename, int species_num)
{
	ifstream file(filename);
	string line;
	getline(file, line);
	while (getline(file, line))
	{
		vector<string> cols;
		split(line, ' ', cols);
		int id=stoi(cols[0]);
		for(int i=0;i<species_num;i++)
		{
			if(!(cols[1+3*i]=="0" && cols[2+3*i]=="0"))
			{
				orfs[id].alt_start_pos[i]=stoi(cols[1+3*i]);
				orfs[id].alt_end_pos[i]=stoi(cols[2+3*i]);
				orfs[id].alt_contig[i]=stoi(cols[3+3*i]);
			}
		}
	}	
}
void print_orf_bounds(vector<Orf>& orfs,string suffix)
{
	ofstream file("orf_bounds"+suffix);
	file << "id";
	for (int i = 0; i < orfs[0].alt_start_pos.size(); i++)
	{
		file << " start" << i << " end" << i << " contig" << i;
	}
	for (int i = 0; i < orfs.size(); i++)
	{
		file << "\n" << i;
		for (int j = 0; j < orfs[i].alt_start_pos.size(); j++)
		{
			file << " " << orfs[i].alt_start_pos[j] << " " << orfs[i].alt_end_pos[j] << " " << orfs[i].alt_contig[j];
		}
	}
}

void get_flank_align(const MultiBlock& block, Orf& orf,int flank_length,int species)
{
	for (int i = 0; i < block.align_seq[0].size(); i++)
	{
		if (block.is_good)
		{
			if (i >= orf.in_start - flank_length && i <= orf.in_end + flank_length)
			{
				if (block.align_seq[1][i] < 0 || block.align_seq[1][i]>6)
				{
					cout << "\nmess: " << i << " " << species << " " << block.align_seq[1][i] << " " << block.is_good;
					getchar();
				}
				orf.focal_flank_align.at(species).push_back(block.align_seq[0][i]);
				orf.flank_align.at(species).push_back(block.align_seq[1][i]);
			}
		}
	}
}

void synteny_match()
{
	int focal_species = 0;
	string suffix = "";
	vector<string> spp = { "Scer","Spar","Smik","Sjur","Skud","Sarb","Suva","Seub" };
	vector<string> genomes = {
		"S288C_reference_sequence_R64-2-1_20150113.fsa" ,
		"Spar.ultrascaf",
		"Smik.ultrascaf",
		"GCA_900290405.1_SacJureiUoM1_genomic.fna",
		"Skud.ultrascaf",
		"GCF_000292725.1_SacArb1.0_genomic.fna",
		"Sbay.ultrascaf",
		"GCF_001298625.1_SEUB3.0_genomic.fna"
	};

	vector<Contig> contigs;
	read_fasta(contigs, genomes[focal_species]);
	//mt19937 rng;
	vector<vector<MultiBlock>> blocks(spp.size());
	for (int i = 0; i < spp.size(); i++)
	{
		read_multiblock(blocks[i], "aligned_blocks" + to_string(i) + "_" + to_string(focal_species) + suffix, 2);
		cout << "\nblocks read: " << blocks[i].size();
		//getchar();
	}
	vector<Orf> orfs;

	//vector<vector<int>> starts = { {1,3,2} };
	get_orfs(orfs, contigs, spp.size(), 7, false);
	remove_duplicate_orfs_by_contig(orfs, contigs);

	map<string, int> id_to_orf;
	for (int i = 0; i < orfs.size(); i++)
	{
		id_to_orf[orfs[i].orf_id] = i;
	}
	for (int i = 0; i < blocks.size(); i++)
	{
		for (int j = 0; j < blocks[i].size(); j++)
		{
			get_block_limits(blocks[i][j]);
		}
		check_blocks_good(blocks[i]);
	}
	cout << "\nget block limits";

	vector<vector<MultiBlock>> rc_blocks = blocks;
	for (int i = 0; i < rc_blocks.size(); i++)
	{
		revc_blocks(rc_blocks[i]);
	}

	vector<vector<double>> sub_rates;
	read_sub_rates(sub_rates, "sub_rates_nongenic");
	for (int i = 0; i < blocks.size(); i++)
	{
		if (i == focal_species)
		{
			continue;
		}
		for (int j = 0; j < blocks[i].size(); j++)
		{
			int orf_index = id_to_orf.at(blocks[i][j].specific_orf);
			orfs[orf_index].frame_seq = vector<int>();
			orfs[orf_index].align_index = vector<int>();
			//cout << "\n" << i << " " << j << " " << blocks[i][j].is_good<<" "<<orf_index<<" "<<blocks[i][j].specific_orf;
			//getchar();
			if (blocks[i][j].is_good)
			{
				if (orfs[orf_index].strand == 0)
				{
					find_orf_in_block(blocks[i][j], orfs[orf_index], 0, i);
					single_compare_orfs(orfs[orf_index], blocks[i][j], i, sub_rates);
					get_flank_align(blocks[i][j], orfs[orf_index], 50, i);
				}
				else
				{
					find_orf_in_block(rc_blocks[i][j], orfs[orf_index], 1, i);
					single_compare_orfs(orfs[orf_index], rc_blocks[i][j], i, sub_rates);
					get_flank_align(blocks[i][j], orfs[orf_index], 50, i);
				}
			}
		}
	}
	cout << "\nfound orf in block";
	print_orf_bounds(orfs,"_synteny");
	cout << "\ncheck good align";
	//getchar();
	//old_check_good_align_(orfs);
	//check_good_align(orfs);
	print_matched_orfs(orfs, "_" + spp[focal_species] + "_synteny_matched" + suffix);
}

void blast_match()
{
	vector<string> spp = { "Scer","Spar","Smik","Sjur","Skud","Sarb","Suva","Seub" };

	vector<Contig> contigs;
	read_fasta(contigs, "S288C_reference_sequence_R64-2-1_20150113.fsa");

	mt19937 rng;
	vector<vector<MultiBlock>> blocks(spp.size());
	for (int i = 1; i < spp.size(); i++)
	{
		read_multiblock(blocks[i], "aligned_blocks_blast_" + to_string(i), 2);
		cout << "\nblocks read: " << blocks[i].size();
	}
	vector<Orf> orfs;
	read_orfs_long(orfs, "orfs_scer_synteny_matched_checked", spp.size());
	get_orf_seqs(orfs, contigs);

	map<string, int> id_to_orf;
	for (int i = 0; i < orfs.size(); i++)
	{
		id_to_orf[orfs[i].orf_id] = i;
	}
	vector<map<string,double>> best_orf_blast_eval(spp.size());
	for (int i = 0; i < blocks.size(); i++)
	{
		for (int j = 0; j < blocks[i].size(); j++)
		{
			get_block_limits(blocks[i][j]);
			if(!best_orf_blast_eval[i].count(blocks[i][j].specific_orf) || blocks[i][j].blast_evalue< best_orf_blast_eval[i].at(blocks[i][j].specific_orf))
			{
				best_orf_blast_eval[i][blocks[i][j].specific_orf] = blocks[i][j].blast_evalue;
			}
		}
		check_blocks_good(blocks[i]);
	}

	vector<vector<MultiBlock>> rc_blocks = blocks;
	for (int i = 0; i < rc_blocks.size(); i++)
	{
		revc_blocks(rc_blocks[i]);
	}

	vector<vector<double>> sub_rates;
	read_sub_rates(sub_rates, "sub_rates_nongenic");
	for (int i = 0; i < blocks.size(); i++)
	{
		for (int j = 0; j < blocks[i].size(); j++)
		{
			int orf_index = id_to_orf.at(blocks[i][j].specific_orf);
			//if(blocks[i][j].specific_orf=="0_334_648_0")
			//{
			//	cout<<"\noi: "<<i<<" "<<j<<" "<<orf_index<<" "<<blocks[i][j].specific_orf<<" "<<blocks[i][j].blast_evalue<<" "<<best_orf_blast_eval[i][blocks[i][j].specific_orf]<<" "<<(blocks[i][j].blast_evalue > best_orf_blast_eval[i][blocks[i][j].specific_orf]);
			//	getchar();
			//}
			if(blocks[i][j].blast_evalue > best_orf_blast_eval[i][blocks[i][j].specific_orf])
			{
				continue;
			}
			orfs[orf_index].frame_seq = vector<int>();
			orfs[orf_index].align_index = vector<int>();
			if (orfs[orf_index].strand == 0)
			{
				Orf dummy = orfs[orf_index];
				find_orf_in_block(blocks[i][j], orfs[orf_index], 0,i);
				if (single_compare_orfs(orfs[orf_index], blocks[i][j], i, sub_rates))
				{
					orfs[orf_index].use_blast_align[i] = 1;
					orfs[orf_index].focal_flank_align[i] = vector<int>();
					orfs[orf_index].flank_align[i] = vector<int>();
					get_flank_align(blocks[i][j], orfs[orf_index], 50, i);
				}
				else
				{
					orfs[orf_index] = dummy;
				}
			}
			else
			{
				Orf dummy = orfs[orf_index];
				find_orf_in_block(rc_blocks[i][j], orfs[orf_index], 1,i);
				if (single_compare_orfs(orfs[orf_index], rc_blocks[i][j], i, sub_rates))
				{
					orfs[orf_index].use_blast_align[i] = 1;
					orfs[orf_index].focal_flank_align[i] = vector<int>();
					orfs[orf_index].flank_align[i] = vector<int>();
					get_flank_align(rc_blocks[i][j], orfs[orf_index], 50, i);
				}
				else
				{
					orfs[orf_index] = dummy;
				}
			}
			//if(blocks[i][j].specific_orf=="0_334_648_0")
			//{
			//	cout<<"\noi: "<<i<<" "<<j<<" "<<orf_index<<" "<<blocks[i][j].specific_orf<<" "<<blocks[i][j].blast_evalue<<" "<<best_orf_blast_eval[i][blocks[i][j].specific_orf]<<" "<<(blocks[i][j].blast_evalue > best_orf_blast_eval[i][blocks[i][j].specific_orf]);
			//	cout<<"\nbi: "<<orfs[orf_index].alt_start_pos[i]<<" "<<orfs[orf_index].alt_end_pos[i];
			//	getchar();
			//}

		}
	}
	print_orf_bounds(orfs,"_blast");
	//check_good_align(orfs);
	print_matched_orfs(orfs, "_scer_blast_matched");

}

void remove_gaps(vector<int> &seq)
{
	vector<int> new_seq;
	for (int i = 0; i < seq.size(); i++)
	{
		if (seq[i] < 4)
		{
			new_seq.push_back(seq[i]);
		}
	}
	seq = new_seq;
}

int count_nongaps(const vector<int>& seq)
{
	int count = 0;
	for (int i = 0; i < seq.size(); i++)
	{
		if (seq[i] != 4)
		{
			count++;
		}
	}
	return count;
}

double water_align(const vector<int> seq0, const vector<int> seq1,double match_bonus,double mismatch_penalty,double gap_penalty)
{
	vector<vector<double>> scores(seq0.size() + 1, vector<double>(seq1.size() + 1));
	double max_score = 0;
	for (int i = 1; i < scores.size(); i++)
	{
		for (int j = 1; j < scores[i].size(); j++)
		{
			double score_match = scores[i - 1][j - 1] + (seq0[i] == seq1[j]) * match_bonus - (seq0[i] != seq1[j]) * mismatch_penalty;
			double score_gap0 = scores[i - 1][j] - gap_penalty;
			double score_gap1 = scores[i][j - 1] - gap_penalty;
			if (score_match > scores[i][j])
			{
				scores[i][j] = score_match;
			}
			if (score_gap0 > scores[i][j])
			{
				scores[i][j] = score_gap0;
			}
			if (score_gap1 > scores[i][j])
			{
				scores[i][j] = score_gap1;
			}
			if (scores[i][j] > max_score)
			{
				max_score = scores[i][j];
			}
		}
	}
	return max_score;
}

void print_fasta(string filename, const vector<int> &seq, string label,bool remove_gaps)
{
	map<int, char> rev_nuc_map = { { 0,'A' },{ 1,'C' },{ 2,'G' },{ 3,'T' },{ 4,'N' },{ 5,'N' },{ 6,'N' } };
	ofstream file(filename);
	file << ">" << label << "\n";
	for (int i = 0; i < seq.size(); i++)
	{
		if (!rev_nuc_map.count(seq[i]))
		{
		}
		else
		{
			if (!remove_gaps || seq[i] != 4)
			{
				file << rev_nuc_map.at(seq[i]);
			}
		}
	}
}

void validate_orf_alignments(vector<Orf> &orfs)
{
	time_t start, finish;
	time(&start);
	for (int i = 0; i < orfs.size(); i++)
	{
		for (int j = 0; j < orfs[i].focal_align_seq.size(); j++)
		{
			int nongaps0 = count_nongaps(orfs[i].focal_align_seq[j]);
			int nongaps1 = count_nongaps(orfs[i].align_seq[j]);
			int align_length = orfs[i].focal_align_seq[j].size();
			//if(orfs[i].use_blast_align[j])
			//{
			//	cout<<"\n"<<i<<" "<<j<<" "<<align_length<<" "<<orfs[i].use_blast_align[j]<<" "<<nongaps0<<" "<<nongaps1<<" "<<orfs[i].align_evalues[j];
			//	getchar();
			//}
			if ( /*orfs[i].use_blast_align[j] &&*/ nongaps0 > 10 &&  nongaps1> 10 && nongaps0 > .4*align_length && nongaps1> .4*align_length && orfs[i].align_evalues[j]>.000001)
			{
				vector<int> seq0 = orfs[i].focal_flank_align[j];// focal_align_seq[j];
				vector<int> seq1 = orfs[i].flank_align[j];// align_seq[j];
				
				remove_gaps(seq0);
				remove_gaps(seq1);
				print_fasta("wwa0.fasta", seq0, "Scer", true);
				print_fasta("wwa1.fasta", seq1, "Scer", true);

				double true_score = water_align(seq0,seq1,5,4,4);
				vector<double> null_scores(1000);
				int better = 0;
				vector<int> null_seq = seq1;
				for (int t = 0; t < null_scores.size(); t++)
				{
					random_shuffle(null_seq.begin(), null_seq.end());
					null_scores[t] = water_align(seq0,null_seq,5,4,4);
					if (true_score > null_scores[t])
					{
						better++;
					}
				}
				orfs[i].align_evalues[j] = 1.0-(double)better / (double)null_scores.size();
				time(&finish);
			}
		}
	}
}

void resolve_intersecting_bounds(vector<Orf> &orfs)
{
	for(int i=0;i<orfs[0].alt_start_pos.size();i++)//iterate through species
	{
		map<int,map<int,vector<int>>> alt_genome_coverage;//chromosome + position to ids of ORFs at that position
		for(int j=0;j<orfs.size();j++)
		{
			if(!(orfs[j].alt_start_pos[i]==0 && orfs[j].alt_end_pos[i]==0))
			{				
				for(int k=orfs[j].alt_start_pos[i];k<=orfs[j].alt_end_pos[i];k++)//iterate through range of homologous region to ORF on comparison species chromosome
				{
					if(alt_genome_coverage[orfs[j].alt_contig[i]].count(k))
					{
						for(int w=0;w<alt_genome_coverage[orfs[j].alt_contig[i]].size();w++)
						{
							int comp_orf_id=alt_genome_coverage[orfs[j].alt_contig[i]][k][w];
							//if(orfs[comp_orf_id].start_pos-orfs[comp_orf_id].start_pos
							cout<<"\nbound inter hit: "<<j<<" "<<comp_orf_id<<" "<<orfs[j].alt_start_pos[i]<<" "<<orfs[j].alt_end_pos[i];
							getchar();
						}
					}
					alt_genome_coverage[orfs[j].alt_contig[i]][k].push_back(j);
				}
			}
		}
	}
}
void validate_blast_matches(const vector<Blast_info> &blasts, vector<Orf> &orfs)
{
	map<string,vector<int>> id_to_blasts;
	for(int i=0;i<blasts.size();i++)
	{
		id_to_blasts[blasts[i].str_id0].push_back(i);
	}
	for(int i=0;i<orfs.size();i++)
	{
		for(int j=0;j<orfs[i].use_blast_align.size();j++)
		{
			string pos_id=orfs[i].orf_id+"_"+to_string(j);
			if(id_to_blasts.count(pos_id))
			{
				double min_eval=1;
				for(int k=0;k<id_to_blasts.at(pos_id).size();k++)
				{
					if(blasts[id_to_blasts.at(pos_id).at(k)].evalue<min_eval)
					{
						min_eval=blasts[id_to_blasts.at(pos_id).at(k)].evalue;
					}
				}
				bool intersects = false;
				for(int k=0;k<id_to_blasts.at(pos_id).size();k++)
				{
					if(blasts[id_to_blasts.at(pos_id).at(k)].evalue<=min_eval)
					{
						int blast_bound0=blasts[id_to_blasts.at(pos_id).at(k)].align1_start;
						int blast_bound1=blasts[id_to_blasts.at(pos_id).at(k)].align1_end;
						if(orfs[i].start_pos<=blast_bound1 && orfs[i].end_pos>=blast_bound0)
						{
							intersects=true;							
						}
					}
				}
				if(intersects)
				{
					orfs[i].homology_validated[j]=true;
				}
			}
		}
	}
	
}
void validate_alignments()
{
	vector<string> spp = { "Scer","Spar","Smik","Sjur","Skud","Sarb","Suva","Seub" };
	vector<Orf> orfs;
	read_orfs_long(orfs, "orfs_scer_blast_matched", spp.size());
	//read_orf_bounds(orfs,"orf_bounds_synteny",spp.size());
	//read_orf_bounds(orfs,"orf_bounds_blast",spp.size());
	//resolve_intersecting_bounds(orfs);
	vector<Blast_info> to_scer_blasts;
	read_blast_results(to_scer_blasts,"blast_matches_to_validate.out",true,true);
	validate_blast_matches(to_scer_blasts,orfs);
	validate_orf_alignments(orfs);
	cout << "\naligned all";
	print_matched_orfs(orfs, "_scer_blast_matched_validated");
}
void print_assigned_homologous_sequence(vector<Orf> &orfs, vector<vector<Contig>> &genomes)
{		
	map<int, string> rev_nuc_map = { { 0,"A" },{ 1,"C" },{ 2,"G" },{ 3,"T" },{ 5,"N" },{ 6,"N" } };

	vector<string> seqs;
	vector<string> ids;
	for(int i=0;i<orfs.size();i++)
	{
		for (int j = 1; j < orfs[i].use_blast_align.size(); j++)
		{
			if(orfs[i].use_blast_align[j])
			{
				seqs.push_back(string());
				ids.push_back(">"+orfs[i].orf_id+"_"+to_string(j));
				for(int k=orfs[i].alt_start_pos[j];k<=orfs[i].alt_end_pos[j];k++)
				{
					seqs.back()+=rev_nuc_map.at(genomes[j][orfs[i].alt_contig[j]].seq[k]);
				}
			}
		}
	}
	ofstream file("blast_matches_to_validate.fasta");
	for(int i=0;i<ids.size();i++)
	{
		file<<ids[i]<<"\n"<<seqs[i]<<"\n";
	}
	string cmd = "blastn -db genome_nuc_Scer -query blast_matches_to_validate.fasta -evalue 10 -word_size 7 -out blast_matches_to_validate.out -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen\"";
	int s = system(cmd.c_str());
}
void resolve_alignments()
{
	vector<string> spp = { "Scer","Spar","Smik","Sjur","Skud","Sarb","Suva","Seub" };
	vector<string> genomes =
	{
		"S288C_reference_sequence_R64-2-1_20150113.fsa" ,
		"Spar.ultrascaf",
		"Smik.ultrascaf",
		"GCA_900290405.1_SacJureiUoM1_genomic.fna",
		"Skud.ultrascaf",
		"GCF_000292725.1_SacArb1.0_genomic.fna",
		"Sbay.ultrascaf",
		"GCF_001298625.1_SEUB3.0_genomic.fna"
	};
	vector<vector<Contig>> contigs(spp.size());
	for (int i = 0; i < spp.size(); i++)
	{
		read_fasta(contigs[i], genomes[i]);
	}
	/////
	/*cout<<"\nseqX:\n"; 
	for(int q=0;q<40;q++)
	{
		cout<<contigs[1][15].seq[590050+q];
	}
	getchar();*/
	////
	vector<Orf> orfs;//orfs_scer_blast_matched
	//read_orfs_long(orfs, "orfs_scer_blast_matched_validated", spp.size());
	read_orfs_long(orfs, "orfs_scer_blast_matched", spp.size());

	//read_orf_bounds(orfs,"orf_bounds_synteny",spp.size());
	read_orf_bounds(orfs,"orf_bounds_blast",spp.size());
	//resolve_intersecting_bounds(orfs);
	print_assigned_homologous_sequence(orfs, contigs);
}


void read_vcfs(vector<Vcf>& vcfs, vector<string>& isolate_names, string filename)
{
	map<char, int> nucmap = { { 'A',0 },{ 'a',0 },{ 'C',1 },{ 'c',1 },{ 'G',2 },{ 'g',2 },{ 'T',3 },{ 't',3 },{ 'N',4 },{ 'n',4 } };

	ifstream file(filename);
	string line;
	while (getline(file, line))
	{
		vector<string> column_data;
		split(line, '\t', column_data);
		if (column_data[0].substr(0, 6) == "#CHROM")
		{
			for (int q = 9; q <= 1019; q++)
			{
				isolate_names.push_back(column_data[q]);
			}
		}
		if (column_data[0].substr(0, 5) == "chrom")
		{
			vcfs.push_back(Vcf());
			vcfs.back().contig = stoi(column_data[0].substr(10)) - 1;
			vcfs.back().pos = stoi(column_data[1]) - 1;
			int entry0 = line.find("AC=") + 3;
			int entry1 = line.find(";AF=");
			vcfs.back().varcount = stoi(line.substr(entry0, entry1 - entry0));
			if (vcfs.back().varcount > 1005)
			{
				vcfs.back().varcount = 2010 - vcfs.back().varcount;
			}
			vcfs.back().variants.push_back(vector<int>());
			for (int k = 0; k < column_data[3].size(); k++)
			{
				vcfs.back().variants.back().push_back(nucmap.at(column_data[3][k]));
			}
			vcfs.back().max_length = vcfs.back().variants.back().size();
			vcfs.back().reference_length = vcfs.back().variants.back().size();
			vector<string> alt_vars;
			split(column_data[4], ',', alt_vars);
			//vcfs.back().variants.push_back(vcfs.back().reference);
			for (int r = 0; r < alt_vars.size(); r++)
			{
				vcfs.back().variants.push_back(vector<int>());
				for (int k = 0; k < alt_vars[r].size(); k++)
				{
					vcfs.back().variants.back().push_back(nucmap.at(alt_vars[r][k]));
				}
				if (vcfs.back().variants.back().size() > vcfs.back().max_length)
				{
					vcfs.back().max_length = vcfs.back().variants.back().size();
				}
			}
			for (int r = 0; r < vcfs.back().variants.size(); r++)
			{
				while (vcfs.back().max_length - vcfs.back().variants[r].size() > 0)
				{
					vcfs.back().variants[r].push_back(4);
				}
			}
			vcfs.back().calls = vector<vector<int>>(1011, vector<int>(2));
			for (int q = 9; q <= 1019; q++)
			{
				string call = column_data[q].substr(0, column_data[q].find(":"));
				vector<string> calls(2);
				calls[0] = call.substr(0, call.find("/"));
				calls[1] = call.substr(call.find("/") + 1, call.size() - 1);
				for (int p = 0; p < 2; p++)
				{
					if (calls[p] == ".")
					{
						vcfs.back().calls[q - 9][p] = -1;
					}
					else
					{
						vcfs.back().calls[q - 9][p] = stoi(calls[p]);
					}
				}
			}
		}
	}
}

void get_major_var(vector<Vcf>& vcfs)
{
	for (int i = 0; i < vcfs.size(); i++)
	{
		vector<int> variant_freqs(vcfs[i].variants.size());
		for (int j = 0; j < vcfs[i].calls.size(); j++)
		{
			for (int k = 0; k < vcfs[i].calls[j].size(); k++)
			{
				if (vcfs[i].calls[j][k] >= 0)
				{
					variant_freqs[vcfs[i].calls[j][k]]++;
				}
			}
		}
		int max_variant_freq = -1;
		for (int j = 0; j < variant_freqs.size(); j++)
		{
			if (variant_freqs[j] > max_variant_freq)
			{
				max_variant_freq = variant_freqs[j];
				vcfs[i].major_var = j;
			}
		}
	}
}

void map_vcfs_to_genome(vector<vector<int>>& vcf_map, const vector<Vcf>& vcfs, const vector<Contig>& contigs)
{
	vcf_map = vector<vector<int>>(contigs.size());
	for (int i = 0; i < contigs.size(); i++)
	{
		vcf_map[i] = vector<int>(contigs[i].seq.size(), -1);
	}
	for (int i = 0; i < vcfs.size(); i++)
	{
		vcf_map[vcfs[i].contig][vcfs[i].pos] = i;
	}
}

void get_genic_status(vector<vector<int>>& genic_status, const vector<Contig>& contigs, const vector<Orf>& genes)
{
	genic_status = vector<vector<int>>(contigs.size());
	for (int i = 0; i < contigs.size(); i++)
	{
		genic_status[i] = vector<int>(contigs[i].seq.size());
	}
	for (int i = 0; i < genes.size(); i++)
	{
		for (int j = genes[i].start_pos; j <= genes[i].end_pos; j++)
		{
			genic_status[genes[i].contig][j] = 1;
		}
	}
}

void get_sub_rates(vector<vector<double>>& sub_rates, const vector<Vcf>& vcfs, const vector<vector<int>>& vcf_map, const vector<vector<int>>& genic_status)
{
	sub_rates = vector<vector<double>>(4, vector<double>(4));
	vector<vector<int>> sub_counts(4, vector<int>(4));
	vector<int> base_counts(4);
	for (int i = 0; i < genic_status.size(); i++)
	{
		for (int j = 0; j < genic_status[i].size(); j++)
		{
			int vcf_id = vcf_map[i][j];
			if (genic_status[i][j] == 0 && vcf_id >= 0 && vcfs[vcf_id].max_length == 1)
			{
				int major_nuc = vcfs[vcf_id].variants[vcfs[vcf_id].major_var][0];
				for (int k = 0; k < vcfs[vcf_id].variants.size(); k++)
				{
					if (vcfs[vcf_id].major_var != k)
					{
						base_counts[major_nuc]++;
						sub_counts[major_nuc][vcfs[vcf_id].variants[k][0]]++;
					}
				}
			}
		}
	}
	for (int i = 0; i < sub_counts.size(); i++)
	{
		for (int j = 0; j < sub_counts[i].size(); j++)
		{
			sub_rates[i][j] = (double)sub_counts[i][j] / (double)base_counts[i];
		}
	}
}


void print_sub_rates(const vector<vector<vector<double>>>& sub_rates, string filename)
{
	ofstream file(filename);
	file << "branch base sub rate";
	for (int i = 0; i < sub_rates.size(); i++)
	{
		for (int j = 0; j < sub_rates[i].size(); j++)
		{
			for(int k=0; k < sub_rates[i][j].size(); k++)
			{
				file << "\n" << i<<" "<<j<<" "<<k<<" "<<sub_rates[i][j][k];				
			}
		}
	}
}

void print_sub_rates(const vector<vector<double>>& sub_rates, string filename)
{
	ofstream file(filename);
	file << "subs 0 1 2 3 4";
	for (int i = 0; i < sub_rates.size(); i++)
	{
		file << "\n" << i;
		for (int j = 0; j < sub_rates[i].size(); j++)
		{
			file << " " << sub_rates[i][j];
		}
	}
}

void get_sub_matrix()
{
	vector<Contig> contigs;
	read_fasta(contigs, "S288C_reference_sequence_R64-2-1_20150113.fsa");
	cout << "\nread contigs: " << contigs.size();

	vector<Vcf> vcfs;
	vector<string> isolates;

	read_vcfs(vcfs, isolates, "/home/acwach/YeastComp/1011Matrix.gvcf");
	cout << "\nvcfs read: " << vcfs.size();

	get_major_var(vcfs);
	cout << "\nget major var";

	vector<vector<int>> vcf_map;
	map_vcfs_to_genome(vcf_map, vcfs, contigs);

	vector<Orf> genes;
	read_annotated_genes(genes, "orf_coding.fasta", true);
	cout << "\nannotated genes read: " << genes.size();

	vector<vector<int>> genic_status;
	get_genic_status(genic_status, contigs, genes);
	cout << "\nget genic status";

	vector<vector<double>> sub_rates;
	get_sub_rates(sub_rates, vcfs, vcf_map, genic_status);
	print_sub_rates(sub_rates, "sub_rates_nongenic");
}

void check_good_align(vector<Orf>& orfs)
{
	//blastn -query seq0 -subject seq1 -word_size 7 -out seq12.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"
	for (int i = 0; i < orfs.size(); i++)
	{
		for (int j = 0; j < orfs[i].focal_align_seq.size(); j++)
		{
			if (count_nongaps(orfs[i].focal_align_seq[j]) > 10 && count_nongaps(orfs[i].align_seq[j]) > 10)
			{
				print_mult_fasta("seq0", { orfs[i].focal_align_seq[j] }, { "Scer" }, false);
				print_mult_fasta("seq1", { orfs[i].align_seq[j] }, { "Alt" }, false);
				string cmd = "blastn -query seq0 -subject seq1 -word_size 7 -out seq12.out -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen\"";
				int s = system(cmd.c_str());
				vector<Blast_info> blasts;
				read_blast_results(blasts, "seq12.out", false, false);
				cout << "\nblasts: " << i << " " << j << " " << orfs[i].is_gene << " " << blasts.size();
				//getchar();
				for (int k = 0; k < blasts.size(); k++)
				{
					if (blasts[k].evalue < orfs[i].align_evalues[j])
					{
						orfs[i].align_evalues[j] = blasts[k].evalue;
					}
				}
			}
		}
	}
}

void synteny_check()
{
	vector<string> spp = { "Scer","Spar","Smik","Sjur","Skud","Sarb","Suva","Seub" };
	vector<Orf> orfs;
	read_matched_orfs(orfs, "orfs_Scer_synteny_matched", spp.size());
	check_good_align(orfs);
	print_matched_orfs(orfs,"_scer_synteny_matched_checked");
	//print_orfs(orfs, "_scer_synteny_matched_checked");
}
void read_repeat_fasta(vector<Contig>& contigs, string filename)
{
	ifstream file(filename);
	string line;
	int contig_id = -1;
	while (getline(file, line))
	{
		if (line.substr(0, 1) == ">")
		{
			contig_id++;
		}
		else
		{
			for (unsigned int i = 0; i < line.size(); i++)
			{
				if (islower(line[i]))
				{
					contigs[contig_id].repeat.push_back(1);
				}
				else
				{
					contigs[contig_id].repeat.push_back(0);
				}
			}
		}
	}
}

void sensu_stricto_blast()
{
	vector<string> spp = { "Scer","Spar","Smik","Sjur","Skud","Sarb","Suva","Seub" };
	vector<string> genomes = {
		"S288C_reference_sequence_R64-2-1_20150113.fsa" ,
		"Spar.ultrascaf",
		"Smik.ultrascaf",
		"GCA_900290405.1_SacJureiUoM1_genomic.fna",
		"Skud.ultrascaf",
		"GCF_000292725.1_SacArb1.0_genomic.fna",
		"Sbay.ultrascaf",
		"GCF_001298625.1_SEUB3.0_genomic.fna"
	};
	vector<Contig> contigs;
	vector<Orf> orfs;
	read_fasta(contigs, "S288C_reference_sequence_R64-2-1_20150113.fsa");
	read_repeat_fasta(contigs, "Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa");
	get_orfs(orfs, contigs, 1, 7, false);
	remove_duplicate_orfs_by_contig(orfs, contigs);
	cout << "\norfs remaining after remove duplicates: " << orfs.size();
	for (int i = 0; i < spp.size(); i++)
	{
		string cmd = "makeblastdb -in " + genomes[i] + " -dbtype nucl -out genome_nuc_" + spp[i];
		int s = system(cmd.c_str());
		cmd = "blastn -db genome_nuc_" + spp[i] + " -query orfs_comp_nuc.fasta -evalue 10 -word_size 7 -out blasts_genome_nuc_" + spp[i] + "_3.out -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen\"";
		s = system(cmd.c_str());
	}
}

void get_mult_align(vector<Orf> &orfs,vector<string> &spp)
{
	map<int, char> rev_nuc_map = { { 0,'A' },{ 1,'C' },{ 2,'G' },{ 3,'T' },{ 4,'-' },{ 5,'N' },{ 6,'N' } };
	ofstream file("orfs_mult_align");
	file << "orf_id";
	map<string, int> labels;
	for (int i = 0; i < spp.size(); i++)
	{
		file << " " << spp[i];
		labels[spp[i]] = i;
	}
	for (int i = 0; i < orfs.size(); i++)
	{
		vector<vector<int>> seqs;
		seqs.push_back(orfs[i].seq);// focal_align_seq[1]);
		seqs.insert(seqs.end(), orfs[i].align_seq.begin()+1, orfs[i].align_seq.end());
		print_mult_fasta("mf0.fas", seqs, spp, false);
		cout << "\norf " << i;

		string cmd = "./muscle3.8.31_i86linux64 -in mf0.fas -fastaout orfs_mult_align.out -quiet";
		int s = system(cmd.c_str());
		vector<string> align_seqs(2);
		vector<vector<int>> mult(spp.size());
		read_muscle(mult,labels, "orfs_mult_align.out");

		file << "\n" << orfs[i].orf_id;
		for (int j = 0; j < mult.size(); j++)
		{
			file << " ";
			for (int k = 0; k < mult[j].size(); k++)
			{
				file << rev_nuc_map.at(mult[j][k]);
			}
		}
	}
}

void mult_align_orthologs()
{
	vector<string> spp = { "Scer","Spar","Smik","Sjur","Skud","Sarb","Suva","Seub" };
	vector<Orf> orfs;		
	read_matched_orfs(orfs, "orfs_scer_blast_matched_validated", spp.size());
	get_mult_align(orfs,spp);
}

void phylo_blast(string blast_type)
{
	vector<string> genome_filenames;
	get_filenames(genome_filenames, "budding_yeast_genomes_filenames.txt",false);
	for (int i = 0; i < genome_filenames.size(); i++)
	{
		//if (genome_filenames[i][0] == '/')
		//{
		//	vector<string> split_filename;
		//	split(genome_filenames[i], '/', split_filename);
			string species_id = genome_filenames[i];//split_filename.back();
			//species_id = species_id.substr(0, species_id.size() - 8);
			//genome_filenames[i].erase(std::remove(genome_filenames[i].begin(), genome_filenames[i].end(), '\n'), genome_filenames[i].end());
			//genome_filenames[i].erase(std::remove(genome_filenames[i].begin(), genome_filenames[i].end(), '\r'), genome_filenames[i].end());
			if (blast_type=="-TBLASTN")
			{
				string msg = "makeblastdb -in budding_yeast_genomes/" + species_id + ".fas -title " + species_id + " -dbtype nucl -out nuc_databases/" + species_id;
				cout << "\n" << msg;
				int s = system(msg.c_str());
				msg = "tblastn -db nuc_databases/" + species_id + " -query orfs_comp.fasta -evalue 1e-1 -out tblastn_align/scer_vs_" + species_id + " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue\"";
				cout << "\n" << msg;
				s = system(msg.c_str());
				//tblastn -db nuc_databases/nakaseomyces_nivariensis -query orfs_comp.fasta -evalue 1e-3 -out tblastn_align/scer_vs_nakaseomyces_nivariensis -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue"
			}
			else if (blast_type=="-TBLASTN_SCRAMBLED")
			{
				//string msg = "tblastn -db nuc_databases/" + species_id + " -query orfs_comp_scrambled.fasta -evalue 1e-1 -out tblastn_align_scrambled/scer_vs_" + species_id + " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue\"";
				string msg = "tblastn -db nuc_databases/" + species_id + " -query orfs_comp_scrambled_and_real.fasta -evalue 1e-1 -out tblastn_align_scrambled_and_real/scer_vs_" + species_id + " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue\"";
				cout << "\n" << msg;
				int s = system(msg.c_str());
			}
			else if (blast_type=="-ORF_BLAST")//yHMPu5000034622_pichia_occidentalis_16_orfs.fas
			{
				string msg = "makeblastdb -in ORFS_332/" + species_id + "_orfs.fas -title " + species_id + " -dbtype prot -out orf_databases/" + species_id;
				cout << "\n" << msg;
				int s = system(msg.c_str());
				msg = "blastp -db orf_databases/" + species_id + " -query orfs_comp.fasta -evalue 1e-3 -out orf_blast_align/scer_vs_" + species_id + " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue\"";
				cout << "\n" << msg;
				s = system(msg.c_str());
			}
			else if(blast_type=="-BLASTP")
			{
				string msg = "makeblastdb -in pep/" + genome_filenames[i] + ".max.pep -title " + species_id + " -dbtype prot -out databases/" + species_id;
				cout << "\n" << msg;
				int s = system(msg.c_str());
				msg = "blastp -db databases/" + species_id + " -query orfs_comp.fasta -evalue 1e-3 -out blastp/scer_vs_" + species_id + " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue\"";
				cout << "\n" << msg;
				s = system(msg.c_str());
			}
		//}
		//candida alibcans special:makeblastdb -in /home/acwach/YeastComp/Genome332/0_332yeast_genomes/Candida_albicans_SC5314_A22_current_default_protein.fasta -title candida_albicans -dbtype prot -out databases/candida_albicans
	}
}

void print_phylo_presence(const vector<vector<double>> &phylo_present, const vector<string> &species_ids,string suffix)
{
	ofstream file("phylo_pres"+suffix);
	file << "orf";
	for (int i = 0; i < species_ids.size(); i++)
	{
		file << " " << species_ids[i];
	}
	for (int i = 0; i < phylo_present.size(); i++)
	{
		file << "\n" << i;
		for (int j = 0; j < phylo_present[i].size(); j++)
		{
			file << " " << phylo_present[i][j];
		}
	}
}

void get_phylo_presence(vector<vector<double>> &phylo_present, vector<vector<BlastMatch>>& blast_match, const map<string,int> &orf_map, const vector<Blast_info> &blasts,double threshold,int genome_id,bool get_opp_info, bool scrambled)
{
	for (int i = 0; i < blasts.size(); i++)
	{
		string orf_id = blasts[i].str_id0;
		if (scrambled)
		{
			orf_id = orf_id.substr(0, orf_id.size() - 2);
		}
		if (orf_map.count(orf_id))
		{
			int orf_index = orf_map.at(orf_id);

			phylo_present[orf_index][genome_id] = blasts[i].evalue;
			blast_match[orf_index][genome_id].opp_gene_id = blasts[i].str_id1.substr(0, blasts[i].str_id1.size() - 2);
			blast_match[orf_index][genome_id].eval = blasts[i].evalue;
			blast_match[orf_index][genome_id].length0 = blasts[i].length0;
			blast_match[orf_index][genome_id].length1 = blasts[i].length1;
			blast_match[orf_index][genome_id].identity = blasts[i].identity;
			blast_match[orf_index][genome_id].mismatch_count = blasts[i].mismatch_count;
			blast_match[orf_index][genome_id].gap_count = blasts[i].gap_count;

			//if (PHYLO_TBLASTN || get_opp_info)
			//{
				blast_match[orf_index][genome_id].contig_id = blasts[i].str_id1;
				blast_match[orf_index][genome_id].opp_start_pos = blasts[i].align1_start;
				blast_match[orf_index][genome_id].opp_stop_pos = blasts[i].align1_end;
			//}
		}
	}
}

void phylo_analyze(string blast_type)
{
	vector<string> genome_filenames;
	get_filenames(genome_filenames, "budding_yeast_genomes_filenames.txt",false);

	vector<Orf> orfs;
	read_orfs(orfs, "orfs_comp", 1);

	map<string, int> orf_map;
	for (int i = 0; i < orfs.size(); i++)
	{
		orf_map[orfs[i].orf_id] = i;
	}
	vector<string> species_ids(genome_filenames.size());
	for (int i = 0; i < genome_filenames.size(); i++)
	{
		species_ids[i]=genome_filenames[i];
	}
	vector<vector<double>> phylo_present(orfs.size(), vector<double>(species_ids.size(), 1));
	vector<vector<double>> phylo_present_scrambled(orfs.size(), vector<double>(species_ids.size(), 1));

	vector<vector<BlastMatch>> blast_match(orfs.size(),vector<BlastMatch>(species_ids.size()));
	for (int i = 0; i < species_ids.size(); i++)
	{
		vector<string> split_filename;
		string blast_filename = "scer_vs_" + species_ids[i];
		if (blast_type=="-TBLASTN")
		{
			blast_filename = "tblastn_align/" + blast_filename;
		}
		else if (blast_type=="-TBLASTN_SCRAMBLED")
		{
			blast_filename = "tblastn_align_scrambled_and_real/" + blast_filename;
		}
		vector<Blast_info> blasts;
		read_blast_results(blasts, blast_filename, false, false);
		cout << "\nblasts read from file " << blast_filename << ":" << blasts.size();
		get_phylo_presence(phylo_present, blast_match, orf_map, blasts, .001, i, false, false);
		get_phylo_presence(phylo_present_scrambled, blast_match, orf_map, blasts, .001, i, false, true);
		string suffix = "";
		if (blast_type=="-TBLASTN")
		{
			suffix = "_tblastn";
		}
		else if (blast_type=="-TBLASTN_SCRAMBLED")
		{
			print_phylo_presence(phylo_present, species_ids, suffix);
			suffix = "_tblastn_scrambled";
			print_phylo_presence(phylo_present_scrambled, species_ids, suffix);
		}
		//if (blast_type=="-TBLASTN_SCRAMBLED")
		//{
		//	print_phylo_presence(phylo_present, species_ids, suffix);
		//}
		/*if (blast_type=="-TBLASTN")
		{
			print_blast_match(blast_match, orf_list, species_ids, "_tblastn");
		}
		else
		{
			print_blast_match(blast_match, orf_list, species_ids, "");
		}*/

		//read gtf file to get coords of best matches
		//read nucleotide fasta file to get actual sequences of matching ORF with flanks
		//for each ORF that overlaps ancient gene in antisense direction, find all ORFs in ancient gene region
		//calculate amino acid similarity with each ORF
		//record stats of best match
	}
}

void phylo_find_orfs() //find ORFs for 332 yeast species by scanning genomes
{
	vector<string> genome_filenames;
	get_filenames(genome_filenames, "budding_yeast_genomes_filenames.txt",false);
	cout << "\ngenomes to read: " << genome_filenames.size();
	for (int i = 0; i < genome_filenames.size(); i++)
	{
		//vector<string> split_filename;
		//split(genome_filenames[i], '/', split_filename);
		string species_id = genome_filenames[i];
		//species_id = species_id.substr(0, species_id.size() - 8);
		cout << "\n" << species_id;
		//getchar();
		//genome_filenames[i].erase(std::remove(genome_filenames[i].begin(), genome_filenames[i].end(), '\n'), genome_filenames[i].end());
		//genome_filenames[i].erase(std::remove(genome_filenames[i].begin(), genome_filenames[i].end(), '\r'), genome_filenames[i].end());
		vector<Contig> contigs;
		vector<Orf> orfs;
		read_fasta(contigs, "budding_yeast_genomes/"+genome_filenames[i]+".fas");
		get_orfs(orfs, contigs, 1, 1, false);
		cout << "\ngenome read: " << genome_filenames[i];
		cout << "\ncontigs read: " << contigs.size();
		cout << "\norfs read: " << orfs.size();

		remove_duplicate_orfs_by_contig(orfs, contigs);
		cout << "\norfs remaining after remove duplicates: " << orfs.size();
		print_orfs_aa_fasta("ORFS_332/"+species_id+"_orfs.fas", orfs);

		//getchar();
	}
}
void print_gene_similar(const vector<GeneSimilar> &gene_similar)
{
	ofstream file("gene_sim");
	file << "id best_match best_evalue best_identity verified dubious pseudogene te all_matches all_evals";
	for (int i = 0; i < gene_similar.size(); i++)
	{
		file <<"\n"<< i << " " << gene_similar[i].best_match << " " << gene_similar[i].evalue << " " << gene_similar[i].identity<<" "<< gene_similar[i].verified_uncharacterized_match << " " << gene_similar[i].dubious_match << " " << gene_similar[i].pseudogene_match << " " << gene_similar[i].te_match<<" ";
		for (int j = 0; j < gene_similar[i].all_matches.size(); j++)
		{
			file<<gene_similar[i].all_matches[j] << ",";
		}
		file << " ";
		for (int j = 0; j < gene_similar[i].all_matches.size(); j++)
		{
			file<<gene_similar[i].all_evalues[j] << ",";
		}
	}
}

void get_gene_similar(vector<GeneSimilar> &gene_similar, const vector<Orf> &orfs, const vector<Blast_info> &blasts, const vector<Orf> &genes)
{
	gene_similar = vector<GeneSimilar>(orfs.size());
	map<string, int> orf_map;
	for (int i = 0; i < orfs.size(); i++)
	{
		orf_map[orfs[i].orf_id] = i;
	}
	map<string, int> gene_map;
	for (int i = 0; i < genes.size(); i++)
	{
		gene_map[genes[i].annotation] = i;
	}
	for (int i = 0; i < blasts.size(); i++)
	{
		int orf_id = orf_map.at(blasts[i].str_id0);
		int gene_id = gene_map.at(blasts[i].str_id1);
		if (!(blasts[i].align0_start == blasts[i].align1_start && blasts[i].align0_end == blasts[i].align1_end && blasts[i].identity==100))
		{
			gene_similar[orf_id].all_matches.push_back(blasts[i].str_id1);
			gene_similar[orf_id].all_evalues.push_back(blasts[i].evalue);
			if (blasts[i].evalue < gene_similar[orf_id].evalue)
			{
				gene_similar[orf_id].evalue = blasts[i].evalue;
				gene_similar[orf_id].best_match = blasts[i].str_id1;
				gene_similar[orf_id].identity = blasts[i].identity;
			}
			if (genes[gene_id].orf_class == "dubious")
			{
				gene_similar[orf_id].dubious_match = true;
			}
			if (genes[gene_id].orf_class == "te")
			{
				gene_similar[orf_id].te_match = true;
			}
			if (genes[gene_id].orf_class == "pseudogene")
			{
				gene_similar[orf_id].pseudogene_match = true;
			}
			if (genes[gene_id].orf_class == "verified" || genes[gene_id].orf_class == "uncharacterized")
			{
				gene_similar[orf_id].verified_uncharacterized_match = true;
			}
		}
	}
}

void self_blast()
{
	vector<Contig> contigs;
	vector<Orf> orfs;
	read_fasta(contigs, "S288C_reference_sequence_R64-2-1_20150113.fsa");
	get_orfs(orfs, contigs, 1, 7, false);
	remove_duplicate_orfs_by_contig(orfs, contigs);
	cout << "\norfs remaining after remove duplicates: " << orfs.size();

	string cmd = "makeblastdb -in orf_genomic_all.fasta -dbtype nucl -out orfs_genomic_all_nuc";
	int s = system(cmd.c_str());

	cmd="blastn -db orfs_genomic_all_nuc -query orfs_comp_nuc.fasta -evalue 1e-3 -word_size 11 -out blasts_scer_orfs_vs_genomic_all -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen\"";
	s=system(cmd.c_str());
	vector<Blast_info> blasts;
	read_blast_results(blasts, "blasts_scer_orfs_vs_genomic_all", false, false);

	vector<Orf> genes;
	read_annotated_genes(genes, "orf_genomic_all.fasta", false);

	vector<GeneSimilar> gene_similar;
	get_gene_similar(gene_similar,orfs, blasts, genes);
	print_gene_similar(gene_similar);
}

void get_hexes(vector<Hex> &hexes)
{
	map<int, char> aa_map = { { 0,'K' },{ 1,'N' },{ 2,'K' },{ 3,'N' },{ 4,'T' },{ 5,'T' },{ 6,'T' },{ 7,'T' },{ 8,'R' },{ 9,'S' },{ 10,'R' },{ 11,'S' },{ 12,'I' },{ 13,'I' },{ 14,'M' },{ 15,'I' },{ 16,'Q' },{ 17,'H' },{ 18,'Q' },{ 19,'H' },{ 20,'P' },{ 21,'P' },{ 22,'P' },{ 23,'P' },{ 24,'R' },{ 25,'R' },{ 26,'R' },{ 27,'R' },{ 28,'L' },{ 29,'L' },{ 30,'L' },{ 31,'L' },{ 32,'E' },{ 33,'D' },{ 34,'E' },{ 35,'D' },{ 36,'A' },{ 37,'A' },{ 38,'A' },{ 39,'A' },{ 40,'G' },{ 41,'G' },{ 42,'G' },{ 43,'G' },{ 44,'V' },{ 45,'V' },{ 46,'V' },{ 47,'V' },{ 48,'X' },{ 49,'Y' },{ 50,'X' },{ 51,'Y' },{ 52,'S' },{ 53,'S' },{ 54,'S' },{ 55,'S' },{ 56,'X' },{ 57,'C' },{ 58,'W' },{ 59,'C' },{ 60,'L' },{ 61,'F' },{ 62,'L' },{ 63,'F' } };
	map<char, int> aa_id = { {'A',0},{'C',1},{'D',2},{'E',3},{'F',4},{'G',5},{'H',6},{'I',7},{'K',8},{'L',9},{'M',10},{'N',11},{'P',12},{'Q',13},{'R',14},{'S',15},{'T',16},{'V',17},{'W',18},{'Y',19},{'X',20} };

	hexes = vector<Hex>(4096);
	for (int i = 0; i < hexes.size(); i++)
	{
		int num = i;
		for (int j = 5; j >= 0; j--)
		{
			hexes[i].seq[j] = num / pow(4, j);
			num -= hexes[i].seq[j] * pow(4, j);
		}
		int aa_trip0 = 0;
		int aa_trip1 = 0;
		for (int k = 0; k < 3; k++)
		{
			hexes[i].trip0 += hexes[i].seq[k] * pow(4, k);
			hexes[i].trip1 += hexes[i].seq[k + 3] * pow(4, k);
			aa_trip0 += hexes[i].seq[2 - k] * pow(4, k);
			aa_trip1 += hexes[i].seq[5 - k] * pow(4, k);

		}
		hexes[i].aa0 = aa_id.at(aa_map.at(aa_trip0));
		hexes[i].aa1 = aa_id.at(aa_map.at(aa_trip1));
		hexes[i].aa0_str = aa_map.at(aa_trip0);
		hexes[i].aa1_str = aa_map.at(aa_trip1);
		hexes[i].diaa = hexes[i].aa0 + hexes[i].aa1 * 20;
	}
}

void get_hex_counts(const vector<Hex> &hexes, vector<int> &hex_counts, vector<double> &diaa_freqs, const vector<Orf> &genes, const set<int> &exclude_triplets, int offset)
{
	vector<int> diaa_counts(diaa_freqs.size());
	int total_count = 0;
	for (int i = 0; i < genes.size(); i++)
	{
		if (!genes[i].splice)
		{
			vector<int> stops(vector<int>(genes[i].seq.size()));
			//id_stops(stops, genes[i].seq);
			for (int j = 3; j < (signed)genes[i].seq.size() - 9 - offset; j += 3)
			{
				bool is_good = true;
				int hex = 0;
				for (int k = 0; k < 6; k++)
				{
					if (!stops[j + k + offset])
					{
						hex += genes[i].seq[j + k + offset] * pow(4, k);
					}
					else
					{
						is_good = false;
					}
				}
				if (is_good && !(exclude_triplets.count(hexes[hex].trip0) || exclude_triplets.count(hexes[hex].trip1)))
				{
					hex_counts[hex]++;
					/*if (hexes[hex].diaa > diaa_counts.size())
					{
						cout << "\n" << hexes[hex].aa0 << " " << hexes[hex].aa1 << " " << hexes[hex].trip0 << " " << hexes[hex].trip1 << " ";
						for (int q = 0; q < hexes[hex].seq.size(); q++)
						{
							cout << hexes[hex].seq[q];
						}
						getchar();
					}
					diaa_counts.at(hexes[hex].diaa)++;*/
					total_count++;
				}
			}
		}
	}
	for (int i = 0; i < hex_counts.size(); i++)
	{
		//hex_freqs[i] = (double)hex_counts[i] / (double)total_count;
	}
	for (int i = 0; i < diaa_counts.size(); i++)
	{
		diaa_freqs[i] = (double)diaa_counts[i] / (double)total_count;
	}
}

void get_hex_freqs(vector<double> &hex_freqs, const vector<int>&hex_counts)
{
	int total = 0;
	for (int i = 0; i < hex_freqs.size(); i++)
	{
		total += hex_counts[i];
	}
	for (int i = 0; i < hex_freqs.size(); i++)
	{
		hex_freqs[i] = (double)hex_counts[i] / (double)total;
	}
}

void get_orf_free_contigs(vector<Contig> &contigs_orfless, const vector<Contig> &contigs, const vector<Orf> &orfs)
{
	vector<Contig> temp = contigs;
	for (int i = 0; i < orfs.size(); i++)
	{
		for (int j = orfs[i].start_pos; j <= orfs[i].end_pos; j++)
		{
			temp[orfs[i].contig].seq[j] = -1;
		}
	}
	for (int i = 0; i < contigs.size(); i++)
	{
		contigs_orfless.push_back(Contig());
		for (int j = 0; j < contigs[i].seq.size(); j++)
		{
			if (temp[i].seq[j] >= 0 && temp[i].seq[j] < 4)
			{
				contigs_orfless.back().seq.push_back(temp[i].seq[j]);
			}
		}
	}
}

void get_nonorfs(vector<Orf> &nonorfs, const vector<Contig> &contigs_orfless, const vector<Orf> &orfs, mt19937 &rng)
{
	for (int i = 0; i < orfs.size(); i++)
	{
		nonorfs.push_back(Orf(1));
		int contig_id = orfs[i].contig;
		int length = orfs[i].end_pos - orfs[i].start_pos;
		uniform_int_distribution<> random_pos(0, contigs_orfless[contig_id].seq.size() - length - 1);
		nonorfs[i].start_pos = random_pos(rng);
		nonorfs[i].end_pos = nonorfs[i].start_pos + length;
		nonorfs[i].strand = orfs[i].strand;
		nonorfs[i].contig = contig_id;
	}
}

void get_hex_scores(vector<Hex> &hexes, const vector<double> &hex_freqs_genes, const vector<double> &hex_freqs_nongenes)
{
	for (int i = 0; i < hexes.size(); i++)
	{
		if (hex_freqs_genes[i] > 0 & hex_freqs_nongenes[i] > 0)
		{
			hexes[i].score = log(hex_freqs_genes[i] / hex_freqs_nongenes[i]);
		}
	}
}

void get_coding_scores(vector<Hex> &hexes, vector<Orf> &orfs, const vector<double> &hex_freq_coding, const vector<double> &hex_freq_noncoding, const vector<double> &diaa_freq_coding, const vector<double> &diaa_freq_noncoding, const set<int> &exclude_triplets, int offset)
{
	for (int i = 0; i < orfs.size(); i++)
	{
		int total_count = 0;
		double coding_score = 0;
		vector<int> stops(orfs[i].seq.size());
		for (int j = 3; j < (signed)orfs[i].seq.size() - 9 - offset; j += 3)
		{
			bool is_good = true;
			int hex = 0;
			for (int k = 0; k < 6; k++)
			{
				if (!stops[j + k + offset])
				{
					hex += orfs[i].seq[j + k + offset] * pow(4, k);
				}
				else
				{
					is_good = false;
				}
			}
			if (is_good && !(exclude_triplets.count(hexes[hex].trip0) || exclude_triplets.count(hexes[hex].trip1)))
			{
				hexes[hex].orf_freq++;
				coding_score += log(hex_freq_coding.at(hex) / hex_freq_noncoding.at(hex));
				total_count++;
			}
		}
		if (total_count > 0)
		{
			orfs[i].coding_score = coding_score / (double)total_count;
		}
	}
}

void print_coding_scores(const vector<Orf> &orfs, string suffix)
{
	ofstream file("coding_scores_" + suffix);
	file << "id score";
	for (int i = 0; i < orfs.size(); i++)
	{
		file << "\n" << i << " " << orfs[i].coding_score;
	}
}

void print_hexes(vector<Hex> &hexes, string suffix)
{
	ofstream file("hex" + suffix);
	file << "index score freq str";
	for (int i = 0; i < hexes[0].pre_count.size(); i++)
	{
		file << " pre" << i << " post" << i;
	}
	for (int i = 0; i < hexes.size(); i++)
	{
		file << "\n" << i << " " << hexes[i].score << " " << hexes[i].orf_freq << " ";
		for (int j = 0; j < hexes[i].seq.size(); j++)
		{
			file << hexes[i].seq[j];
		}
		for (int j = 0; j < hexes[i].pre_count.size(); j++)
		{
			file << " " << hexes[i].pre_count[j] << " " << hexes[i].post_count[j];
		}
	}
}

void coding_scores()
{
	mt19937 rng;

	vector<Contig> contigs;
	read_fasta(contigs, "S288C_reference_sequence_R64-2-1_20150113.fsa");

	vector<Orf> genes;
	read_annotated_genes(genes, "orf_coding.fasta", true);
	cout<<"\ngenes read: "<<genes.size();

	vector<Orf> orfs;
	read_orfs(orfs, "orfs_comp", 8);
	get_orf_seqs(orfs, contigs);
	cout << "\norfs read: " << orfs.size();

	set<int> exclude_triplets;
	exclude_triplets.insert({ 48,50,56 });

	vector<Hex> hexes;
	get_hexes(hexes);

	vector<int> hex_counts(4096);
	vector<double> hex_freqs(4096);
	vector<double> diaa_freqs(400);
	get_hex_counts(hexes, hex_counts, diaa_freqs, genes, exclude_triplets, 0);
	get_hex_freqs(hex_freqs, hex_counts);

	vector<Contig> contigs_orfless;
	get_orf_free_contigs(contigs_orfless, contigs, genes);
	cout<<"\ncontigs read: "<<contigs_orfless.size();

	vector<Orf> nonorfs;
	get_nonorfs(nonorfs, contigs_orfless, genes, rng);
	sort(nonorfs.begin(), nonorfs.end());
	get_orf_seqs(nonorfs, contigs_orfless);
	cout << "\nnonorfs count:" << nonorfs.size();

	vector<int> control_hex_counts(4096);
	vector<double> control_hex_freqs(4096);
	vector<double> control_diaa_freqs(4096);

	get_hex_counts(hexes, control_hex_counts, control_diaa_freqs, nonorfs, exclude_triplets, 0);
	get_hex_freqs(control_hex_freqs,control_hex_counts);

	cout << "\nget coding scores";
	get_hex_scores(hexes, hex_freqs, control_hex_freqs);

	get_coding_scores(hexes, orfs, hex_freqs, control_hex_freqs, diaa_freqs, control_diaa_freqs, exclude_triplets, 0);

	print_coding_scores(orfs, "");
	print_hexes(hexes, "");
}

void print_pop_align(const vector<vector<vector<int>>> &pop_align, const vector<int> absent, int suffix, const vector<string> &isolates)
{
	ofstream file("ORF_COMP_POP/pop_align" + to_string(suffix));
	file << "pop_align" << suffix;
	for (int i = 0; i < pop_align.size(); i++)
	{
		file << "\n" << isolates[i];
		if (!absent[i])
		{
			for (int k = 0; k < 2; k++)
			{
				file << " ";
				for (int j = 0; j < pop_align[i][k].size(); j++)
				{
					file << pop_align[i][k][j];
				}
			}
		}
	}
}


void match_proto_vcf(vector<vector<int>> &align_positions,vector<Orf> &protogenes, const vector<Contig> &contigs, const vector<Vcf> &vcfs, const vector<vector<int>> vcf_map, const vector<string> &isolates)
{
	for (int i = 0; i < protogenes.size(); i++)
	{
		vector<vector<vector<int>>> mult_align(1011, vector<vector<int>>(2));
		vector<int> absent(1011);
		for (int j = protogenes[i].start_pos; j <= protogenes[i].end_pos; j++)
		{
			int vcf_id = vcf_map.at(protogenes[i].contig).at(j);
			if (vcf_id >= 0)
			{
				for (int k = 0; k < vcfs.at(vcf_id).variants[0].size() && k+j<=protogenes[i].end_pos; k++)
				{
					if (vcfs.at(vcf_id).variants[0][k] != 4)
					{
						align_positions[i].push_back(j-protogenes[i].start_pos+k);
					}
					else
					{
						align_positions[i].push_back(-1);
					}
				}
				for (int v = 0; v < vcfs.at(vcf_id).calls.size(); v++)
				{
					if (!absent[v])
					{
						for (int p = 0; p < vcfs[vcf_id].calls[v].size(); p++)
						{
							int call = vcfs.at(vcf_id).calls.at(v).at(p);
							if (call == -1)//inconclusive call
							{
								absent.at(v) = 1;
							}
							else
							{
								for (int k = 0; k < vcfs.at(vcf_id).variants[call].size() && k+j<=protogenes[i].end_pos; k++)
								{
									mult_align[v][p].push_back(vcfs.at(vcf_id).variants[call][k]);
								}
								//mult_align[v][p].insert(mult_align[v][p].end(), vcfs.at(vcf_id).variants[call].begin(), vcfs.at(vcf_id).variants[call].end());
							}
						}
					}
				}
				j += vcfs.at(vcf_id).reference_length - 1;
			}
			else
			{
				for (int v = 0; v < mult_align.size(); v++)
				{
					if (!absent[v])
					{
						for (int p = 0; p < mult_align[v].size(); p++)
						{
							mult_align[v][p].push_back(contigs[protogenes[i].contig].seq[j]);
						}
					}
				}
				align_positions[i].push_back(j - protogenes[i].start_pos);
			}
		}
		print_pop_align(mult_align, absent, i, isolates);
	}
}
void print_align_positions(const vector<vector<int>> &align_positions)
{
	ofstream file("align_positions");
	for (int i = 0; i < align_positions.size(); i++)
	{
		file << "\n" << i;
		for (int j = 0; j < align_positions[i].size(); j++)
		{
			file << " " << align_positions[i][j];
		}
	}
}

void strain_setup()
{
	vector<Contig> contigs;
	read_fasta(contigs, "S288C_reference_sequence_R64-2-1_20150113.fsa");

	vector<Orf> orfs;
	read_orfs(orfs, "orfs_comp", 8);
	cout<<"\norfs read: "<<orfs.size();
	getchar();
	
	vector<Vcf> vcfs;
	vector<string> isolates;

	read_vcfs(vcfs, isolates, "1011Matrix.gvcf");
	cout << "\nvcfs read: " << vcfs.size();
	getchar();
	
	vector<vector<int>> vcf_map;
	map_vcfs_to_genome(vcf_map, vcfs, contigs);

	//get_pop_align_to_ref(orfs, vcfs, vcf_map, contigs);
	vector<vector<int>> align_positions(orfs.size());
	match_proto_vcf(align_positions, orfs, contigs, vcfs, vcf_map, isolates);
	print_align_positions(align_positions);

}

void read_mult_align(vector<vector<string>> &mult_align, string filename)
{
	ifstream file(filename);
	string line;
	getline(file,line);
	while (getline(file, line))
	{
		vector<string> columns;
		split(line, ' ', columns);
		mult_align.push_back(vector<string>());
		for(int i=1;i<columns.size();i++)
		{
			mult_align.back().push_back(columns[i]);
		}
	}
}

void get_matched_alignment(vector<string> &match_align, vector<string> &mult_align, vector<int> &included_species, vector<vector<int>> &match_frame)
{
	int codon_pos=0;
	match_align=vector<string>(1+included_species.size());
	for(int i=0;i<mult_align[0].size()-2;i++)
	{
		if(mult_align[0][i]!='-')
		{
			if(codon_pos==0)
			{
				bool all_good=true;
				for(int j=0;j<match_frame.size();j++)
				{
					//cout<<"\nmaz: "<<mult_align[0].size()-2<<" "<<match_frame[j].size();
					//getchar();
					if(match_frame[j].at(i)==1 && match_frame[j].at(i+1)==1 && match_frame[j].at(i+2)==1 && mult_align[0][i]!='-' && mult_align[0][i+1]!='-' && mult_align[0][i+2]!='-')
					{
						
					}
					else
					{
						all_good=false;
					}
				}
				if(all_good)
				{
					match_align[0]+=mult_align[0][i];
					match_align[0]+=mult_align[0][i+1];
					match_align[0]+=mult_align[0][i+2];					
					for(int j=0;j<included_species.size();j++)
					{
						match_align[j+1]+=mult_align[included_species[j]][i];
						match_align[j+1]+=mult_align[included_species[j]][i+1];
						match_align[j+1]+=mult_align[included_species[j]][i+2];					
					}
				}
			}
			codon_pos++;
			if(codon_pos==3)
			{
				codon_pos=0;
			}
		}
	}
	//cout<<"\nicm: "<<included_species.size()<<" "<<match_align.size();
	//getchar();
	//ofstream file("matched_alignment");
	//file<<">spp0\n"<<match_align[0];
	//for(int i=1;i<match_align.size();i++)
	//{
	//	file<<"\n>spp"<<included_species[i-1]<<"\n"<<match_align[i];
	//}
}

bool seq_has_n(vector<int> &seq)
{
	for(int i=0;i<seq.size();i++)
	{
		if(seq[i]==5)
		{
			return true;
		}
	}
	return false;
}

int count_match_frames(vector<int> match_frame)
{
	int count=0;
	for(int i=0;i<match_frame.size();i++)
	{
		if(match_frame[i]==1)
		{
			count++;
		}
	}
	return count;
}

void make_clean_alignment_files()
{
	vector<string> spp = { "Scer","Spar","Smik","Sjur","Skud","Sarb","Suva","Seub" };
	ofstream file("clean_alignments");
	
	file<<"orf_id";
	for(int i=0;i<spp.size();i++)
	{
		file<<" "<<spp[i];
	}
	// load sequences and reading frame for each ORF/species
	vector<Orf> orfs;
	read_orfs_long(orfs, "orfs_scer_blast_matched_validated", spp.size());
	cout << "\norfs read: " << orfs.size()<<" "<<orfs[0].orf_id;
	vector<vector<string>> mult_align;
	read_mult_align(mult_align, "orfs_mult_align");
	cout<<"\nmult align read: "<<mult_align.size();
	//getchar();
	for(int i=0;i<orfs.size();i++)
	{
		vector<int> included_species;
		for(int j=1;j<orfs[i].align_evalues.size();j++)
		{
			if(orfs[i].align_evalues[j]<.01 && (orfs[i].use_blast_align[j]==0 || orfs[i].homology_validated[j]==1) && orfs[i].match_frame[j].size()>0 &&!seq_has_n(orfs[i].align_seq[j]) && (double)count_match_frames(orfs[i].match_frame[j])/(double)orfs[i].match_frame[j].size() >.4 )
			{
				included_species.push_back(j);
			}
			/*if(orfs[i].orf_id=="4_565602_565799_0")
			{
				cout<<"\norf: "<<i<<" "<<orfs[i].align_evalues[j]<<" "<<orfs[i].use_blast_align[j]<<" "<<orfs[i].homology_validated[j]<<" "<<orfs[i].match_frame[j].size()<<" "<<included_species.size()<<" "<<seq_has_n(orfs[i].align_seq[j]);
				getchar();
			}*/

		}
		if(included_species.size()==0)
		{
			continue;
		}
		cout<<"\nincluded species: "<<included_species.size()<<" "<<i<<" "<<orfs[i].orf_id;
		
		/*if(orfs[i].orf_id=="4_565602_565799_0")
		{
			cout<<"\norf: "<<i;//<<" "<<spp_included;
			getchar();
		}*/


		//getchar();
		vector<vector<int>> match_frame(included_species.size());
		for(int j=0;j<included_species.size();j++)
		{
			int spp=included_species[j];
			vector<int> match_frame_by_nuc;
			//cout<<"\nmfbn: "<<spp<<" "<<match_frame_by_nuc.size()<<" "<<orfs[i].match_frame[spp].size()<<" "<<orfs[i].align_seq[spp].size();
			//cout<<"\nB "<<orfs[i].match_frame[spp].size()<<" "<<spp<<"\n";
			//getchar();
			for(int k=0;k<orfs[i].match_frame[spp].size();k++)
			{
				//cout<<orfs[i].align_seq[spp][k];
				if(orfs[i].align_seq[spp][k]<4)
				{
					match_frame_by_nuc.push_back(orfs[i].match_frame[spp][k]);
				}
			}
			//cout<<"\nmfbn: "<<spp<<" "<<match_frame_by_nuc.size()<<" "<<orfs[i].match_frame[spp].size();
			//getchar();
			//cout<<"\nC "<<match_frame_by_nuc.size();
			//getchar();
			int nuc_pos=0;
			for(int k=0;k<mult_align[i][spp].size();k++)
			{
				if(mult_align[i][spp][k]!='-')
				{
					match_frame[j].push_back(match_frame_by_nuc.at(nuc_pos));
					nuc_pos++;
				}
				else
				{
					match_frame[j].push_back(0);
				}
			}
		}
		vector<string> matched_alignment;
		//cout<<"A";
		//getchar();
		get_matched_alignment(matched_alignment, mult_align[i],included_species,match_frame);
		vector<string> matched_alignment_allspp(spp.size());
		matched_alignment_allspp[0]=matched_alignment[0];
		for(int s=0;s<included_species.size();s++)
		{
			matched_alignment_allspp[included_species[s]]=matched_alignment[s+1];
		}
		file<<"\n"<<orfs[i].orf_id;
		for(int s=0;s<spp.size();s++)
		{
			file<<" "<<matched_alignment_allspp[s];
		}
		//cout<<"\nprint matched alignment for ORF "<<i<<" "<<orfs[i].orf_id;
		//getchar();
	}
	


	// do mult alignment of all orthologous regions
	// id codons with shared reading frame for all orthologs
	// make alignment of those codons
	// do ancestral reconstruction
	// get sub rates on each branch irrespective of codon position
	// calc dnds

}

void read_prank_alignments(vector<string> &prank_alignments, string filename)
{
	ifstream file(filename);
	string line;
	while (getline(file, line))
	{
		if(line[0]!='>')
		{
			prank_alignments.push_back(line);
		}
		else
		{
			cout<<"\n"<<line;
		}
	}
}
void read_clean_alignments(vector<vector<string>> &clean_alignments, vector<string> &orf_ids, string filename)
{
	ifstream file(filename);
	string line;
	getline(file,line);
	while (getline(file, line))
	{
		vector<string> columns;
		split(line, ' ', columns);
		clean_alignments.push_back(vector<string>());
		orf_ids.push_back(columns[0]);
		for(int i=1;i<columns.size();i++)
		{
			clean_alignments.back().push_back(columns[i]);
		}
	}
}

void get_prank_subs(vector<vector<vector<int>>> &subs,vector<vector<int>> &bases,vector<string> seqs)
{
	map<char, int> nucmap = {{ 'A',0 },{ 'C',1 },{ 'G',2 },{ 'T',3 }};
	vector<int> branch_ancestors =   {1,1,5,5,13,13,3,3,7,7,9, 9,11,11};
	vector<int> branch_descendants = {0,2,4,6,12,14,1,5,3,8,7,10, 9,13};
	for(int i=0;i<seqs[0].size();i++)
	{
		for(int j=0;j<branch_ancestors.size();j++)
		{
			bases.at(j).at(nucmap.at(seqs[branch_ancestors[j]][i]))++;
			subs.at(j).at(nucmap.at(seqs[branch_ancestors[j]][i]))[nucmap.at(seqs[branch_descendants[j]][i])]++;
		}
	}
}

void make_joined_clean_alignment()
{
	vector<string> spp = { "Scer","Spar","Smik","Sjur","Skud","Sarb","Suva","Seub" };
	vector<Orf> orfs;
	read_orfs(orfs, "orfs_comp", 1);
	map<string,int> orf_id_map;
	for(int i=0;i<orfs.size();i++)
	{
		orf_id_map[orfs[i].orf_id]=i;
	}
	vector<vector<string>> clean_alignments;
	vector<string> orf_ids;
	read_clean_alignments(clean_alignments,orf_ids,"clean_alignments");
	cout<<"\nread alignments: "<<clean_alignments.size();
	int full_set=0;
	vector<vector<vector<int>>> subs(14, vector<vector<int>>(4,vector<int>(4)));
	vector<vector<int>> bases(14,vector<int>(4));
	vector<string> joined_alignments(spp.size());
	for(int i=0;i<clean_alignments.size();i++)
	{
		int spp_included=0;
		ofstream file("clean_orf.fasta");
		for(int j=0;j<clean_alignments[i].size();j++) 
		{
			if(clean_alignments[i][j].size()>0)
			{
				spp_included++;
			}
		}
		if(spp_included==8 && orfs.at(orf_id_map.at(orf_ids[i])).is_gene=="X" && orfs.at(orf_id_map.at(orf_ids[i])).overlaps_gene=="X")
		{
			full_set++;
			for(int k=0;k<8;k++)
			{
				joined_alignments[k]+=clean_alignments[i][k];
			}
		}
	}
	ofstream file("joined_alignment");
	for(int k=0;k<joined_alignments.size();k++)
	{
		file<<">"<<spp[k]<<"\n"<<joined_alignments[k]<<"\n";
	}
	
	string cmd = "prank -d=joined_alignment -showanc -keep -t=joined_alignment_phyml_tree.txt -once";
	int s = system(cmd.c_str());
	cmd="mv output.anc.dnd joined_tree.nwk";
	s=system(cmd.c_str());
}

void tree_dnds(DNDS_reporter &dnds,const vector<string> &prank_alignments,const vector<vector<double>> &sub_rates, const vector<int> branch_ancestors, const vector<int> branch_descendants)
{
	
	map<char, int> nucmap = {{ 'A',0 },{ 'C',1 },{ 'G',2 },{ 'T',3 }};
	//vector<int> branch_ancestors =   {1,1,5,5,13,13,3,3,7,7,9, 9,11,11};
	//vector<int> branch_descendants = {0,2,4,6,12,14,1,5,3,8,7,10, 9,13};
	map<int, char> aa_map = { { 0,'K' },{ 1,'N' },{ 2,'K' },{ 3,'N' },{ 4,'T' },{ 5,'T' },{ 6,'T' },{ 7,'T' },{ 8,'R' },{ 9,'S' },{ 10,'R' },{ 11,'S' },{ 12,'I' },{ 13,'I' },{ 14,'M' },{ 15,'I' },{ 16,'Q' },{ 17,'H' },{ 18,'Q' },{ 19,'H' },{ 20,'P' },{ 21,'P' },{ 22,'P' },{ 23,'P' },{ 24,'R' },{ 25,'R' },{ 26,'R' },{ 27,'R' },{ 28,'L' },{ 29,'L' },{ 30,'L' },{ 31,'L' },{ 32,'E' },{ 33,'D' },{ 34,'E' },{ 35,'D' },{ 36,'A' },{ 37,'A' },{ 38,'A' },{ 39,'A' },{ 40,'G' },{ 41,'G' },{ 42,'G' },{ 43,'G' },{ 44,'V' },{ 45,'V' },{ 46,'V' },{ 47,'V' },{ 48,'X' },{ 49,'Y' },{ 50,'X' },{ 51,'Y' },{ 52,'S' },{ 53,'S' },{ 54,'S' },{ 55,'S' },{ 56,'X' },{ 57,'C' },{ 58,'W' },{ 59,'C' },{ 60,'L' },{ 61,'F' },{ 62,'L' },{ 63,'F' } };

	vector<int> codon_pos;
	vector<int> seq_pos;
	int current_codon_pos = 0;
	int current_seq_pos = 0;
	for (int j = 0; j < prank_alignments[0].size(); j++)
	{
		codon_pos.push_back(current_codon_pos);
		current_codon_pos++;
		if (current_codon_pos == 3)
		{
			current_codon_pos = 0;
		}
		seq_pos.push_back(current_seq_pos);
		current_seq_pos++;
	}
	for(int i=0;i<prank_alignments[0].size();i++)
	{
		for(int j=0;j<branch_ancestors.size();j++)
		{
			int h = codon_pos[i];
			vector<int> cod0 = { nucmap.at(prank_alignments[branch_ancestors[j]][i - h]), nucmap.at(prank_alignments[branch_ancestors[j]][i + 1 - h]) , nucmap.at(prank_alignments[branch_ancestors[j]][i + 2 - h]) };
			vector<int> cod1 = { nucmap.at(prank_alignments[branch_descendants[j]][i - h]), nucmap.at(prank_alignments[branch_descendants[j]][i + 1 - h]) , nucmap.at(prank_alignments[branch_descendants[j]][i + 2 - h]) };
			char aa0 = aa_map.at(cod0[0] * 4 * 4 + cod0[1] * 4 + cod0[2]);
			for (int k = 0; k < 4; k++)
			{
				if (cod0[h] != k)
				{
					vector<int> pos_cod = cod0;
					pos_cod[h] = k;
					char pos_aa = aa_map.at(pos_cod[0] * 4 * 4 + pos_cod[1] * 4 + pos_cod[2]);

					if (pos_aa == 'X' || pos_aa == 'M')
					{
						continue;
					}
					
					if (pos_aa == aa0)
					{
						dnds.syn_base[cod0[h] * 4 + k]++;
						dnds.syn_exp+=sub_rates[cod0[h]][k];
						if (cod1[h] == k)
						{
							dnds.syn_obs++;
							dnds.syn_subs[cod0[h] * 4 + k]++;
						}
					}
					else
					{
						dnds.nonsyn_base[cod0[h] * 4 + k]++;
						dnds.nonsyn_exp+=sub_rates[cod0[h]][k];
						if (cod1[h] == k)
						{
							dnds.nonsyn_obs++;
							dnds.nonsyn_subs[cod0[h] * 4 + k]++;
						}
					}
				}
			}

		}
	}
}

void print_dnds(const vector<DNDS_reporter> &dnds, string filename)
{
	ofstream file(filename);
	file << "orf_id pass syn_exp nonsyn_exp syn_obs nonsyn_obs";
	for (int i = 0; i < dnds[0].nonsyn_base.size(); i++)
	{
		file << " nonsyn_base" << i;
	}
	for (int i = 0; i < dnds[0].nonsyn_subs.size(); i++)
	{
		file << " nonsyn_subs" << i;
	}
	for (int i = 0; i < dnds[0].syn_base.size(); i++)
	{
		file << " syn_base" << i;
	}
	for (int i = 0; i < dnds[0].syn_subs.size(); i++)
	{
		file << " syn_subs" << i;
	}
	for (int i = 0; i < dnds.size(); i++)
	{
		file << "\n" << dnds[i].orf_id << " " << dnds[i].pass << " " << dnds[i].syn_exp << " " << dnds[i].nonsyn_exp << " " << dnds[i].syn_obs << " " << dnds[i].nonsyn_obs;
		for (int j = 0; j < dnds[i].nonsyn_base.size(); j++)
		{
			file << " " << dnds[i].nonsyn_base[j];
		}
		for (int j = 0; j < dnds[i].nonsyn_subs.size(); j++)
		{
			file << " " << dnds[i].nonsyn_subs[j];
		}
		for (int j = 0; j < dnds[i].syn_base.size(); j++)
		{
			file << " " << dnds[i].syn_base[j];
		}
		for (int j = 0; j < dnds[i].syn_subs.size(); j++)
		{
			file << " " << dnds[i].syn_subs[j];
		}
	}
}


void individual_dnds()
{
	vector<string> spp = { "Scer","Spar","Smik","Sjur","Skud","Sarb","Suva","Seub" };
	vector<Orf> orfs;
	read_orfs(orfs, "orfs_comp", 1);
	map<string,int> orf_id_map;
	for(int i=0;i<orfs.size();i++)
	{
		orf_id_map[orfs[i].orf_id]=i;
	}
	
	vector<vector<string>> clean_alignments;
	vector<string> orf_ids;
	read_clean_alignments(clean_alignments,orf_ids,"clean_alignments");
	cout<<"\nread alignments: "<<clean_alignments.size()<<" "<<orf_ids.size();
	

	int full_set=0;
	vector<vector<vector<int>>> subs(14, vector<vector<int>>(4,vector<int>(4)));
	vector<vector<int>> bases(14,vector<int>(4));

	vector<vector<double>> sub_rates;
	read_sub_rates(sub_rates, "sub_rates_nongenic");

	//vector<string> prank_alignments;
	//read_prank_alignments(prank_alignments,"output.anc.fas");
	//cout<<"\n"<<prank_alignments.size();
	//getchar();
	//get_prank_subs(subs,bases,prank_alignments);
	vector<DNDS_reporter> dnds_reporter(clean_alignments.size());
	for(int i=0;i<clean_alignments.size();i++)
	{
		int spp_included=0;
		ofstream file("clean_orf.fasta");
		for(int j=0;j<clean_alignments[i].size();j++) 
		{
			if(clean_alignments[i][j].size()>0)
			{
				spp_included++;
				file<<">"<<spp[j]<<"\n"<<clean_alignments[i][j]<<"\n";
			}
		}
		file.close();
		/*if(orf_ids[i]=="4_565602_565799_0")
		{
			cout<<"\norf: "<<i<<" "<<spp_included;
			getchar();
		}*/
		vector<int> branch_ancestors;
		vector<int> branch_descendants;
		if(spp_included==8) //has homologs in all species
		{
			branch_ancestors =   {1,1,5,5,13,13,3,3,7,7,9, 9,11,11};
			branch_descendants = {0,2,4,6,12,14,1,5,3,8,7,10, 9,13};
			full_set++;
		}
		else if(spp_included==7 && clean_alignments[i][1].size()==0) //missing paradoxus
		{
			branch_ancestors =   {1,3,3,5,11,5,7,1,7,9,11,9};
			branch_descendants = {0,2,4,6,12,1,5,3,8,7,10,11};
			full_set++;			
		}
		else if(spp_included==7 && clean_alignments[i][4].size()==0) //kud
		{
			branch_ancestors =   {1,3,1,7,5,3,5,9,7,11,9 ,11};
			branch_descendants = {0,1,2,3,4,5,6,7,8,10,11,12};
			full_set++;			
		}
		else if(spp_included==7 && clean_alignments[i][2].size()==0) //mik
		{
			branch_ancestors =   {1,3,1,5,3,7,5,9,7,9 };
			branch_descendants = {0,1,2,3,4,5,6,8,9,10};
			full_set++;			
		}
		else if(spp_included==6 && clean_alignments[i][2].size()==0 && clean_alignments[i][3].size()==0) //missing mikatei/jurei
		{
			branch_ancestors = 	 {1,3,1,5,3,7,5,7};
			branch_descendants = {0,1,2,3,4,6,7,8};
		}
		else if(spp_included==6 && clean_alignments[i][3].size()==0 && clean_alignments[i][5].size()==0)//missing jur/arb
		{
			branch_ancestors =	 {1,3,1,5,3,7,5,9,7,9};
			branch_descendants = {0,1,2,3,4,5,6,8,9,10};
		}
		else if(spp_included==6 && clean_alignments[i][3].size()==0 && clean_alignments[i][4].size()==0)//missing jur/kud
		{
			branch_ancestors =	 {1,3,1,5,3,7,5,9,7,9};
			branch_descendants = {0,1,2,3,4,5,6,8,9,10};
		}
		else if(spp_included==5 && clean_alignments[i][3].size()==0 && clean_alignments[i][4].size()==0 && clean_alignments[i][5].size()==0) //missing jur kud arb
		{
			branch_ancestors   = {1,3,1,5,3,6,5,7};
			branch_descendants = {0,1,2,3,4,6,7,8};
		}
		if(branch_ancestors.size()>0 && orfs.at(orf_id_map.at(orf_ids[i])).overlaps_gene=="X")
		{
			string cmd="rm output.anc.fas";
			int s = system(cmd.c_str());

			cmd = "prank -d=clean_orf.fasta -t=joined_tree.nwk -showanc -showevents -once -prunetree -keep";
			s = system(cmd.c_str());
			
			//if(spp_included==7)
			//{
			//	cout<<"got8";
			//	getchar();
			//}

			vector<string> prank_alignments;
			if(file_exists("output.anc.fas"))
			{
				read_prank_alignments(prank_alignments,"output.anc.fas");
				dnds_reporter[i].orf_id=orf_ids[i];
				tree_dnds(dnds_reporter[i],prank_alignments,sub_rates, branch_ancestors,branch_descendants);
			}			
		}
	}
	//cout<<"A";
	//getchar();
	vector<DNDS_reporter> dnds_reporter_by_orf(orfs.size());
	for(int i=0;i<dnds_reporter_by_orf.size();i++)
	{
		dnds_reporter_by_orf[i].orf_id = orfs[i].orf_id;
	}
	//cout<<"B: "<<orf_ids.size()<<" "<<dnds_reporter.size()<<" "<<dnds_reporter_by_orf.size();
	//getchar();

	for(int i=0;i<dnds_reporter.size();i++)
	{
		/*if(!orf_id_map.count(orf_ids[i]))
		{
			cout<<"\nxyz: "<<i<<" "<<orf_ids[i];
			getchar();
		}*/
		dnds_reporter_by_orf[orf_id_map.at(orf_ids[i])]=dnds_reporter[i];
	}
	//cout<<"C";
	//getchar();

	print_dnds(dnds_reporter_by_orf,"tree_dnds");

	//print_dnds(dnds_reporter,"tree_dnds");
	/*vector<vector<vector<double>>> sub_rates(14, vector<vector<double>>(4,vector<double>(4)));
	for(int i=0;i<sub_rates.size();i++)
	{
		for(int j=0;j<sub_rates[i].size();j++)
		{
			for(int k=0;k<sub_rates[i][j].size();k++)
			{
				sub_rates[i][j][k]=(double)subs[i][j][k]/(double)bases[i][j];
			}
		}
	}
	
	cout<<"\nfull set: "<<full_set;
	print_sub_rates(sub_rates,"branch_sub_rates");*/
	//make_clean_alignment_files();	
	//cmd = "/home/acwach/YeastComp/prank/bin/prank -o=blocks/output_ancestors" + suffix + to_string(i) + " -d=blocks/w" + suffix + to_string(i) + ".fasta -t=blocks/wm" + suffix + to_string(i) + ".out_phyml_tree.txt -showanc -showevents -once";
	//s = system(cmd.c_str());
}

void read_align_positions(vector<vector<int>> &align_positions, string filename)
{
	ifstream file(filename);
	string line;
	getline(file, line);
	while (getline(file, line))
	{
		vector<string> column_data;
		split(line, ' ', column_data);
		int index = stoi(column_data[0]);
		for (int i = 1; i < column_data.size(); i++)
		{
			align_positions[index].push_back(stoi(column_data[i]));
		}
	}
}

void read_pop_align(vector<vector<vector<int>>> &pop_align, int suffix)
{
	ifstream file("/home/acwach/Synteny/ORF_COMP_POP/pop_align" + to_string(suffix));
	string line;
	getline(file, line);
	int i = 0;
	while (getline(file, line))
	{
		vector<string> column_data;
		split(line, ' ', column_data);
		if (column_data.size() > 1)
		{
			for (int k = 0; k < 2; k++)
			{
				for (int j = 0; j < column_data[k + 1].size(); j++)
				{
					pop_align[i][k].push_back(column_data[k + 1][j] - '0');
				}
			}
		}
		i++;
	}
}

void get_matched_sequences(vector<vector<vector<int>>> &matched_sequences,const Orf &orf, const vector<vector<vector<int>>> &pop_align, const vector<int> &align_positions)
{
	matched_sequences = vector<vector<vector<int>>>(pop_align.size(), vector<vector<int>>(2));
	for (int i = 0; i < pop_align.size(); i++)
	{
		for (int j = 0; j < pop_align[i].size(); j++)
		{
			for (int k = 0; k < pop_align[i][j].size(); k++)
			{
				if (align_positions[k] >= 0)
				{
					matched_sequences[i][j].push_back(pop_align[i][j][k]);
				}
			}
			if (orf.strand == 1)
			{
				reverse_complement(matched_sequences[i][j]);
			}
		}
	}
}

void get_variant_freqs(vector<vector<int>> &variant_freqs, const vector<vector<vector<int>>> &matched_sequences)
{
	for (int i = 0; i < matched_sequences.size() && variant_freqs.size()==0; i++)
	{
		for (int j = 0; j < matched_sequences[i].size(); j++)
		{
			if (matched_sequences[i][j].size() > 0)
			{
				variant_freqs = vector<vector<int>>(matched_sequences[i][j].size(), vector<int>(5));
				//cout << "vfs: " << variant_freqs.size();
				//getchar();
				break;
			}
		}
	}
	for (int i = 0; i < matched_sequences.size(); i++)
	{
		for (int j = 0; j < matched_sequences[i].size(); j++)
		{
			for (int k = 0; k < matched_sequences[i][j].size(); k++)
			{
				int nuc = matched_sequences[i][j][k];
				if (k >= variant_freqs.size())
				{
					cout << "\nbad k: " << k << " " << variant_freqs.size();
					getchar();
				}
				if (nuc<0 || nuc>=variant_freqs[k].size())
				{
					cout << "\nbad nuc: " << nuc;
					getchar();
				}
				variant_freqs[k][nuc]++;
			}
		}
	}
}

void get_consensus(vector<int> &consensus, const vector<vector<int>> &variant_freqs)
{
	consensus= vector<int>(variant_freqs.size(), -1);
	for (int i = 0; i < variant_freqs.size(); i++)
	{
		int major_nuc_freq = -1;
		for (int j = 0; j < variant_freqs[i].size(); j++)
		{
			if (variant_freqs[i][j] > major_nuc_freq)
			{
				major_nuc_freq = variant_freqs[i][j];
				consensus[i] = j;
			}
		}
	}
}

void get_syn_nonsyn_vars(Orf &orf, const vector<vector<int>> &variant_freqs,const vector<vector<double>> &sub_rates)
{
	map<int, char> aa_map = { { 0,'K' },{ 1,'N' },{ 2,'K' },{ 3,'N' },{ 4,'T' },{ 5,'T' },{ 6,'T' },{ 7,'T' },{ 8,'R' },{ 9,'S' },{ 10,'R' },{ 11,'S' },{ 12,'I' },{ 13,'I' },{ 14,'M' },{ 15,'I' },{ 16,'Q' },{ 17,'H' },{ 18,'Q' },{ 19,'H' },{ 20,'P' },{ 21,'P' },{ 22,'P' },{ 23,'P' },{ 24,'R' },{ 25,'R' },{ 26,'R' },{ 27,'R' },{ 28,'L' },{ 29,'L' },{ 30,'L' },{ 31,'L' },{ 32,'E' },{ 33,'D' },{ 34,'E' },{ 35,'D' },{ 36,'A' },{ 37,'A' },{ 38,'A' },{ 39,'A' },{ 40,'G' },{ 41,'G' },{ 42,'G' },{ 43,'G' },{ 44,'V' },{ 45,'V' },{ 46,'V' },{ 47,'V' },{ 48,'X' },{ 49,'Y' },{ 50,'X' },{ 51,'Y' },{ 52,'S' },{ 53,'S' },{ 54,'S' },{ 55,'S' },{ 56,'X' },{ 57,'C' },{ 58,'W' },{ 59,'C' },{ 60,'L' },{ 61,'F' },{ 62,'L' },{ 63,'F' } };
	vector<int> consensus;
	get_consensus(consensus, variant_freqs);
	int codon_pos = 0;
	double total_mut_exp = 0;
	int total_muts = 0;
	for (int i = 0; i < variant_freqs.size(); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (j != consensus[i])
			{
				int old_codon_hash = -1;
				int new_codon_hash = -1;
				if (codon_pos == 0)
				{
					if (consensus[i] >= 4 || consensus[i + 1] >= 4 || consensus[i + 2] >= 4)
					{
						continue;
					}
					old_codon_hash = consensus[i] * 4 * 4 + consensus[i + 1] * 4 + consensus[i + 2];
					new_codon_hash = j * 4 * 4 + consensus[i + 1] * 4 + consensus[i + 2];
				}
				if (codon_pos == 1)
				{
					if (consensus[i-1] >= 4 || consensus[i] >= 4 || consensus[i + 1] >= 4)
					{
						continue;
					}
					old_codon_hash = consensus[i - 1] * 4 * 4 + consensus[i] * 4 + consensus[i + 1];
					new_codon_hash= consensus[i - 1] * 4 * 4 + j * 4 + consensus[i + 1];
				}
				if (codon_pos == 2)
				{
					if (consensus[i-2] >= 4 || consensus[i - 1] >= 4 || consensus[i] >= 4)
					{
						continue;
					}
					old_codon_hash = consensus[i - 2] * 4 * 4 + consensus[i - 1] * 4 + consensus[i];
					new_codon_hash = consensus[i - 2] * 4 * 4 + consensus[i - 1] * 4 + j;
				}
				if (!aa_map.count(old_codon_hash) || !aa_map.count(new_codon_hash))
				{
					cout << "\nhash mismatch: " << old_codon_hash << " " << new_codon_hash<<" "<<codon_pos<<" "<<consensus[i]<<consensus[i+1]<<consensus[i+2]<<" "<<j<<" "<<i<<" "<<consensus.size()<<" "<<codon_pos<<" "<<orf.orf_id;
					getchar();
				}
				if (aa_map.at(old_codon_hash) == aa_map.at(new_codon_hash))
				{
					orf.syn_exp += sub_rates[consensus[i]][j];
					total_mut_exp+= sub_rates[consensus[i]][j];
					if (variant_freqs[i][j] > 0)
					{
						orf.syn_vars++;
						total_muts++;
					}
				}
				else
				{
					orf.nonsyn_exp+=sub_rates[consensus[i]][j];
					total_mut_exp+= sub_rates[consensus[i]][j];
					if (variant_freqs[i][j] > 0)
					{
						orf.nonsyn_vars++;
						total_muts++;
					}
				}
			}
		}
		codon_pos++;
		if (codon_pos == 3)
		{
			codon_pos = 0;
		}
	}
	if (total_mut_exp == 0)
	{
		orf.syn_exp = 0;
		orf.nonsyn_exp = 0;
	}
	else
	{
		orf.syn_exp *= (double)total_muts / total_mut_exp;
		orf.nonsyn_exp *= (double)total_muts / total_mut_exp;
	}
}

void get_anti_syn_nonsyn_vars(vector<DNDS_reporter2> &dnds, Orf& orf, const vector<vector<int>>& variant_freqs, const vector<vector<double>>& sub_rates, int rel_frame)
{
	map<int, string> rev_nuc_map = { { 0,"A" },{ 1,"C" },{ 2,"G" },{ 3,"T" } };
	map<int, char> aa_map = { { 0,'K' },{ 1,'N' },{ 2,'K' },{ 3,'N' },{ 4,'T' },{ 5,'T' },{ 6,'T' },{ 7,'T' },{ 8,'R' },{ 9,'S' },{ 10,'R' },{ 11,'S' },{ 12,'I' },{ 13,'I' },{ 14,'M' },{ 15,'I' },{ 16,'Q' },{ 17,'H' },{ 18,'Q' },{ 19,'H' },{ 20,'P' },{ 21,'P' },{ 22,'P' },{ 23,'P' },{ 24,'R' },{ 25,'R' },{ 26,'R' },{ 27,'R' },{ 28,'L' },{ 29,'L' },{ 30,'L' },{ 31,'L' },{ 32,'E' },{ 33,'D' },{ 34,'E' },{ 35,'D' },{ 36,'A' },{ 37,'A' },{ 38,'A' },{ 39,'A' },{ 40,'G' },{ 41,'G' },{ 42,'G' },{ 43,'G' },{ 44,'V' },{ 45,'V' },{ 46,'V' },{ 47,'V' },{ 48,'X' },{ 49,'Y' },{ 50,'X' },{ 51,'Y' },{ 52,'S' },{ 53,'S' },{ 54,'S' },{ 55,'S' },{ 56,'X' },{ 57,'C' },{ 58,'W' },{ 59,'C' },{ 60,'L' },{ 61,'F' },{ 62,'L' },{ 63,'F' } };
	vector<int> rc_nuc = { 3,2,1,0 };
	vector<int> consensus;
	get_consensus(consensus, variant_freqs);
	vector<int> codon_pos(variant_freqs.size());
	int cur_codon_pos = 0;
	for (int i = 0; i < codon_pos.size(); i++)
	{
		codon_pos[i] = cur_codon_pos;
		cur_codon_pos++;
		if (cur_codon_pos == 3)
		{
			cur_codon_pos = 0;
		}
	}
	double total_mut_exp = 0;
	int total_muts = 0;
	vector<int> anti_codon_pos;// = { 0,1,2 };

	if (rel_frame == 3)
	{
		anti_codon_pos = { 0,1,2 };
	}
	if (rel_frame == 2)
	{
		anti_codon_pos = { 1,2,0 };
	}
	if (rel_frame == 1)
	{
		anti_codon_pos = { 2,0,1 };
	}

	for (int i = 5; i < (signed)variant_freqs.size()-5; i++)
	{
		if (consensus[i - codon_pos[i]] >= 4 || consensus[i + 1 - codon_pos[i]] >= 4 || consensus[i + 2 - codon_pos[i]] >= 4||
			consensus[i - anti_codon_pos[codon_pos[i]]]>=4|| consensus[i + 1 - anti_codon_pos[codon_pos[i]]]>=4|| consensus[i + 2 - anti_codon_pos[codon_pos[i]]]>=4||
			consensus[i-1]>=4||consensus[i]>=4||consensus[i+1]>=4)
		{
			continue;
		}

		vector<int> cod0 = { consensus[i - codon_pos[i]],consensus[i + 1 - codon_pos[i]] ,consensus[i + 2 - codon_pos[i]] };
		vector<int> anti_cod0 = { consensus[i-anti_codon_pos[codon_pos[i]]],   consensus[i + 1 - anti_codon_pos[codon_pos[i]]], consensus[i + 2 - anti_codon_pos[codon_pos[i]]] };

		char aa0 = aa_map.at(cod0[0] * 4 * 4 + cod0[1] * 4 + cod0[2]);
		char anti_aa0 = aa_map.at(rc_nuc.at(anti_cod0[2]) * 4 * 4 + rc_nuc.at(anti_cod0[1]) * 4 + rc_nuc.at(anti_cod0[0]));

		for (int k = 0; k < 4; k++)
		{
			if (k != consensus[i])
			{
				vector<int> pos_cod = cod0;
				pos_cod[codon_pos[i]] = k;
				char pos_aa = aa_map.at(pos_cod[0] * 4 * 4 + pos_cod[1] * 4 + pos_cod[2]);

				vector<int> pos_anti_cod = anti_cod0;
				pos_anti_cod[anti_codon_pos[codon_pos[i]]] = k;
				char pos_anti_aa = aa_map.at(rc_nuc.at(pos_anti_cod[2]) * 4 * 4 + rc_nuc.at(pos_anti_cod[1]) * 4 + rc_nuc.at(pos_anti_cod[0]));

				if (pos_anti_aa != anti_aa0)
				{
					continue;
				}

				if (pos_aa == 'X' || pos_anti_aa == 'X' || pos_aa == 'M')
				{
					continue;
				}
				dnds.push_back(DNDS_reporter2());

				dnds.back().context = rev_nuc_map.at(consensus[i - 1]) + rev_nuc_map.at(consensus[i]) + rev_nuc_map.at(consensus[i + 1]);
				dnds.back().base = rev_nuc_map.at(consensus[i]);
				dnds.back().sub = rev_nuc_map.at(k);
				dnds.back().orf_id = orf.orf_id;
				dnds.back().codon = rev_nuc_map.at(cod0[0]) + rev_nuc_map.at(cod0[1]) + rev_nuc_map.at(cod0[2]);
				dnds.back().codon_pos = codon_pos[i];
				dnds.back().anticodon = rev_nuc_map.at(rc_nuc.at(anti_cod0[2])) + rev_nuc_map.at(rc_nuc.at(anti_cod0[1])) + rev_nuc_map.at(rc_nuc.at(anti_cod0[0]));

				if (pos_aa == aa0) 
				{
					dnds.back().syn = 1;
				}
				if (variant_freqs[i][k] > 0)
				{
					dnds.back().actual = 1;
				}
			}
		}
	}

}

void assess_pop_align(Orf &protogene, vector<vector<vector<int>>> &pop_align, int strand)
{
	protogene.has_atg = 0;
	protogene.absent_count = 0;
	protogene.no_atg = 0;
	for (int i = 0; i < pop_align.size(); i++)
	{
		if (pop_align[i][0].size() == 0)
		{
			protogene.absent_count += 2;
		}
		else
		{
			for (int j = 0; j < pop_align[i].size(); j++)
			{
				vector<int> nuc;
				vector<int> seq = pop_align[i][j];
				remove_gaps(seq);
				bool has_atg = false;
				if (strand == 1)
				{
					reverse_complement(seq);
				}
				if (seq[0] == 0 && seq[1] == 3 && seq[2] == 2)
				{
					protogene.has_atg++;
					has_atg = true;
				}
				else
				{
					protogene.no_atg++;
				}
				bool has_stop = false;
				bool intermediate_stop = false;
				for (int k = 0; k < seq.size(); k += 3)
				{
					if (seq[k] == 3 && (
						(seq[k + 1] == 0 && seq[k + 2] == 0) ||
						(seq[k + 1] == 2 && seq[k + 2] == 0) ||
						(seq[k + 1] == 0 && seq[k + 2] == 2)
						))
					{
						if (k > seq.size() - 10)
						{
							has_stop = true;
						}
						else
						{
							intermediate_stop = true;
						}
					}
				}
				protogene.has_stop += has_stop;
				protogene.inter_stop += intermediate_stop;
				if (has_atg && has_stop && !intermediate_stop)
				{
					protogene.conserved_orf++;
				}
			}
		}
	}
}

void print_pop_info(const vector<Orf> &orfs)
{
	ofstream file("pop_info");
	file << "index orf_id has_atg absent_count no_atg has_stop inter_stop conserved_orf";
	for (int i = 0; i < orfs.size();i++)
	{
		file << "\n" << i << " " << orfs[i].orf_id << " " << orfs[i].has_atg << " " << orfs[i].absent_count << " " << orfs[i].no_atg << " " << orfs[i].has_stop << " " << orfs[i].inter_stop << " " << orfs[i].conserved_orf;
	}
}

void print_dnds2(const vector<DNDS_reporter2>& dnds, string filename)
{
	ofstream file(filename);
	file << "orf_id context base sub syn actual codon anticodon codon_pos";
	for (int i = 0; i < dnds.size(); i++)
	{
		file << "\n" << dnds[i].orf_id << " " << dnds[i].context << " " << dnds[i].base << " " << dnds[i].sub << " " << dnds[i].syn << " " << dnds[i].actual<<" "<<dnds[i].codon<<" "<<dnds[i].anticodon<<" "<<dnds[i].codon_pos;
	}
}

void strain_analyze()
{
	vector<vector<double>> sub_rates;
	read_sub_rates(sub_rates, "sub_rates_nongenic");

	vector<Contig> contigs;
	read_fasta(contigs, "S288C_reference_sequence_R64-2-1_20150113.fsa");

	vector<Orf> orfs;
	read_orfs(orfs, "orfs_comp", 1);

	vector<vector<int>> align_positions(orfs.size());
	read_align_positions(align_positions, "align_positions");
	vector<vector<vector<int>>> codon_trans(orfs.size());
	cout << "\nread align positions";
	vector<vector<DNDS_reporter2>> dnds(3);
	for (int i = 0; i < orfs.size(); i++)
	{
		vector<vector<vector<int>>> pop_align(1011, vector<vector<int>>(2));
		vector<vector<int>> variant_freqs;
		read_pop_align(pop_align, i);
		vector<vector<vector<int>>> matched_sequences;
		get_matched_sequences(matched_sequences, orfs[i], pop_align, align_positions[i]);

		if (matched_sequences[0][0].size() != 0 && matched_sequences[0][0].size() != orfs[i].end_pos - orfs[i].start_pos + 1)
		{
			cout << "\n" << i << " " << orfs[i].end_pos - orfs[i].start_pos + 1 << " " << matched_sequences[0][0].size() << " ";
			getchar();
		}
		get_variant_freqs(variant_freqs, matched_sequences);
		vector<int> consensus;
		get_consensus(consensus, variant_freqs);
		//get_codon_trans(codon_trans[i], matched_sequences,consensus );
		get_syn_nonsyn_vars(orfs[i], variant_freqs, sub_rates);
		for (int k = 0; k < 3; k++)
		{
			get_anti_syn_nonsyn_vars(dnds[k], orfs[i], variant_freqs, sub_rates, k+1);
		}
		assess_pop_align(orfs[i], pop_align, orfs[i].strand);
	}

	print_pop_info(orfs);
	print_matched_orfs(orfs, "_pop");
	for (int k = 0; k < 3; k++)
	{
		print_dnds2(dnds[k], "anti_pnps_pop_"+to_string(k+1));
	}	
}

/*bool check_good_dnds(const Orf& orf)
{
	if (orf.match_frame[1].size() == 0)
	{
		return false;
	}
	else if (orf.match_frame[1].size() != orf.focal_align_seq[1].size())
	{
		cout<<"ERROR: "<<orf.orf_id<<" "<<orf.match_frame[1].size()<<" "<<orf.focal_align_seq[1].size();
		getchar();
		return false;
	}
	return true;
}*/


void do_dnds(/*vector<DNDS_reporter>& dnds, */vector<DNDS_reporter2>& dnds2, const vector<Orf>& orfs)
{
	map<int, string> rev_nuc_map = { { 0,"A" },{ 1,"C" },{ 2,"G" },{ 3,"T" } };
	map<int, char> aa_map = { { 0,'K' },{ 1,'N' },{ 2,'K' },{ 3,'N' },{ 4,'T' },{ 5,'T' },{ 6,'T' },{ 7,'T' },{ 8,'R' },{ 9,'S' },{ 10,'R' },{ 11,'S' },{ 12,'I' },{ 13,'I' },{ 14,'M' },{ 15,'I' },{ 16,'Q' },{ 17,'H' },{ 18,'Q' },{ 19,'H' },{ 20,'P' },{ 21,'P' },{ 22,'P' },{ 23,'P' },{ 24,'R' },{ 25,'R' },{ 26,'R' },{ 27,'R' },{ 28,'L' },{ 29,'L' },{ 30,'L' },{ 31,'L' },{ 32,'E' },{ 33,'D' },{ 34,'E' },{ 35,'D' },{ 36,'A' },{ 37,'A' },{ 38,'A' },{ 39,'A' },{ 40,'G' },{ 41,'G' },{ 42,'G' },{ 43,'G' },{ 44,'V' },{ 45,'V' },{ 46,'V' },{ 47,'V' },{ 48,'X' },{ 49,'Y' },{ 50,'X' },{ 51,'Y' },{ 52,'S' },{ 53,'S' },{ 54,'S' },{ 55,'S' },{ 56,'X' },{ 57,'C' },{ 58,'W' },{ 59,'C' },{ 60,'L' },{ 61,'F' },{ 62,'L' },{ 63,'F' } };
	for (int i = 0; i < orfs.size(); i++)
	{
		//dnds[i].orf_id = orfs[i].orf_id;
		if (orfs[i].match_frame[1].size()==0)//!check_good_dnds(orfs[i]))
		{
			continue;
		}
		vector<int> codon_pos;
		vector<int> seq_pos;
		int current_codon_pos = 0;
		int current_seq_pos = 0;
		for (int j = 0; j < orfs[i].focal_align_seq[1].size(); j++)
		{
			if (orfs[i].focal_align_seq[1][j] < 4)
			{
				codon_pos.push_back(current_codon_pos);
				current_codon_pos++;
				if (current_codon_pos == 3)
				{
					current_codon_pos = 0;
				}
				seq_pos.push_back(current_seq_pos);
				current_seq_pos++;
			}
			else
			{
				seq_pos.push_back(-1);
				codon_pos.push_back(4);
			}
		}
		for (int j = 0; j < (signed)orfs[i].focal_align_seq[1].size() - 2; j++)
		{
			if (seq_pos[j] > 4 && seq_pos[j] < current_seq_pos - 5)
			{
				int h = codon_pos[j];
				vector<int> cod0 = { orfs[i].focal_align_seq[1][j - codon_pos[j]],orfs[i].focal_align_seq[1][j + 1 - codon_pos[j]] ,orfs[i].focal_align_seq[1][j + 2 - codon_pos[j]] };
				vector<int> cod1 = { orfs[i].align_seq[1][j - codon_pos[j]],orfs[i].align_seq[1][j + 1 - codon_pos[j]], orfs[i].align_seq[1][j + 2 - codon_pos[j]] };
				if (orfs[i].focal_align_seq[1][j] > 3 || orfs[i].focal_align_seq[1][j - 1] > 3 || orfs[i].focal_align_seq[1][j + 1] > 3)
				{
					continue;
				}
				if (cod0[0] > 3 || cod0[1] > 3 || cod0[2] > 3 || cod1[0] > 3 || cod1[1] > 3 || cod1[2] > 3)
				{
					continue;
				}
				if (orfs[i].match_frame[1][j - codon_pos[j]] != 1 || orfs[i].match_frame[1][j + 1 - codon_pos[j]] != 1 || orfs[i].match_frame[1][j + 2 - codon_pos[j]] != 1)
				{
					continue;
				}
				int sum_diff = (cod0[0] != cod1[0]) + (cod0[1] != cod1[1]) + (cod0[2] != cod1[2]);
				if (sum_diff > 1)
				{
					continue;
				}
				if (!aa_map.count(cod0[0] * 4 * 4 + cod0[1] * 4 + cod0[2]))
				{
					cout << "\nerror: " << cod0[0] << " " << cod0[1] << " " << cod0[2];
					getchar();
				}
				char aa0 = aa_map.at(cod0[0] * 4 * 4 + cod0[1] * 4 + cod0[2]);
				for (int k = 0; k < 4; k++)
				{
					if (cod0[h] != k)
					{
						vector<int> pos_cod = cod0;
						pos_cod[h] = k;
						char pos_aa = aa_map.at(pos_cod[0] * 4 * 4 + pos_cod[1] * 4 + pos_cod[2]);

						if (pos_aa == 'X' || pos_aa == 'M')
						{
							continue;
						}
						dnds2.push_back(DNDS_reporter2());
						dnds2.back().context = rev_nuc_map.at(orfs[i].focal_align_seq[1][j - 1]) + rev_nuc_map.at(orfs[i].focal_align_seq[1][j]) + rev_nuc_map.at(orfs[i].focal_align_seq[1][j + 1]);
						dnds2.back().base = rev_nuc_map.at(orfs[i].focal_align_seq[1][j]);
						dnds2.back().sub = rev_nuc_map.at(k);
						dnds2.back().orf_id = orfs[i].orf_id;
						if (pos_aa == aa0)
						{
							//dnds[i].syn_base[cod0[h] * 4 + k]++;
							dnds2.back().syn = 1;
							if (cod1[h] == k)
							{
								//dnds[i].syn_obs++;
								//dnds[i].syn_subs[cod0[h] * 4 + k]++;
								dnds2.back().actual = 1;
							}
						}
						else
						{
							dnds2.back().syn = 0;
							//dnds[i].nonsyn_base[cod0[h] * 4 + k]++;
							if (cod1[h] == k)
							{
								//dnds[i].nonsyn_obs++;
								//dnds[i].nonsyn_subs[cod0[h] * 4 + k]++;
								dnds2.back().actual = 1;
							}
						}
					}
				}

			}
		}
	}
}

void do_anti_dnds(/*vector<DNDS_reporter>& dnds, */vector<DNDS_reporter2>& dnds2, const vector<Orf>& orfs, int rel_frame, bool anti_synonymous)
{
	map<int, string> rev_nuc_map = { { 0,"A" },{ 1,"C" },{ 2,"G" },{ 3,"T" } };
	map<int, char> aa_map = { { 0,'K' },{ 1,'N' },{ 2,'K' },{ 3,'N' },{ 4,'T' },{ 5,'T' },{ 6,'T' },{ 7,'T' },{ 8,'R' },{ 9,'S' },{ 10,'R' },{ 11,'S' },{ 12,'I' },{ 13,'I' },{ 14,'M' },{ 15,'I' },{ 16,'Q' },{ 17,'H' },{ 18,'Q' },{ 19,'H' },{ 20,'P' },{ 21,'P' },{ 22,'P' },{ 23,'P' },{ 24,'R' },{ 25,'R' },{ 26,'R' },{ 27,'R' },{ 28,'L' },{ 29,'L' },{ 30,'L' },{ 31,'L' },{ 32,'E' },{ 33,'D' },{ 34,'E' },{ 35,'D' },{ 36,'A' },{ 37,'A' },{ 38,'A' },{ 39,'A' },{ 40,'G' },{ 41,'G' },{ 42,'G' },{ 43,'G' },{ 44,'V' },{ 45,'V' },{ 46,'V' },{ 47,'V' },{ 48,'X' },{ 49,'Y' },{ 50,'X' },{ 51,'Y' },{ 52,'S' },{ 53,'S' },{ 54,'S' },{ 55,'S' },{ 56,'X' },{ 57,'C' },{ 58,'W' },{ 59,'C' },{ 60,'L' },{ 61,'F' },{ 62,'L' },{ 63,'F' } };
	vector<int> rc_nuc = { 3,2,1,0 };
	for (int i = 0; i < orfs.size(); i++)
	{
		//dnds[i].orf_id = orfs[i].orf_id;
		if (orfs[i].match_frame[1].size()==0)
		{
			continue;
		}
		vector<int> codon_pos;
		vector<int> seq_pos;
		int current_codon_pos = 0;
		int current_seq_pos = 0;
		for (int j = 0; j < orfs[i].focal_align_seq[1].size(); j++)
		{
			if (orfs[i].focal_align_seq[1][j] < 4)
			{
				codon_pos.push_back(current_codon_pos);
				current_codon_pos++;
				if (current_codon_pos == 3)
				{
					current_codon_pos = 0;
				}
				seq_pos.push_back(current_seq_pos);
				current_seq_pos++;
			}
			else
			{
				seq_pos.push_back(-1);
				codon_pos.push_back(4);
			}
		}
		for (int j = 0; j < (signed)orfs[i].focal_align_seq[1].size() - 2; j++)
		{
			vector<int> anti_codon_pos;
			if (rel_frame == 3)
			{
				anti_codon_pos = { 0,1,2 };
			}
			if (rel_frame == 2)
			{
				anti_codon_pos = { 1,2,0 };
			}
			if (rel_frame == 1)
			{
				anti_codon_pos = { 2,0,1 };
			}
			if (seq_pos[j] > 4 && seq_pos[j] < current_seq_pos - 5 /*&& relevant_positions_match*/)
			{
				int h = codon_pos[j];
				vector<int> cod0 = { orfs[i].focal_align_seq[1][j - codon_pos[j]],orfs[i].focal_align_seq[1][j + 1 - codon_pos[j]] ,orfs[i].focal_align_seq[1][j + 2 - codon_pos[j]] };
				vector<int> cod1 = { orfs[i].align_seq[1][j - codon_pos[j]],orfs[i].align_seq[1][j + 1 - codon_pos[j]], orfs[i].align_seq[1][j + 2 - codon_pos[j]] };
				vector<int> anti_cod0 = { orfs[i].focal_align_seq[1][j - anti_codon_pos[codon_pos[j]]],orfs[i].focal_align_seq[1][j + 1 - anti_codon_pos[codon_pos[j]]] ,orfs[i].focal_align_seq[1][j + 2 - anti_codon_pos[codon_pos[j]]] };
				vector<int> anti_cod1 = { orfs[i].align_seq[1][j - anti_codon_pos[codon_pos[j]]],orfs[i].align_seq[1][j + 1 - anti_codon_pos[codon_pos[j]]], orfs[i].align_seq[1][j + 2 - anti_codon_pos[codon_pos[j]]] };
				if (orfs[i].focal_align_seq[1][j] > 3|| orfs[i].focal_align_seq[1][j-1] > 3|| orfs[i].focal_align_seq[1][j + 1] > 3)
				{
					continue;
				}

				if (cod0[0] > 3 || cod0[1] > 3 || cod0[2] > 3 || cod1[0] > 3 || cod1[1] > 3 || cod1[2] > 3 ||
					anti_cod0[0] > 3 || anti_cod0[1] > 3 || anti_cod0[2] > 3 || anti_cod1[0] > 3 || anti_cod1[1] > 3 || anti_cod1[2] > 3)
				{
					continue;
				}
				if (orfs[i].match_frame[1][j - codon_pos[j]]!=1||orfs[i].match_frame[1][j + 1 - codon_pos[j]]!=1 || orfs[i].match_frame[1][j+2-codon_pos[j]]!=1)
				{
					continue;
				}
				if (orfs[i].match_frame[1][j - anti_codon_pos[codon_pos[j]]] != 1 || orfs[i].match_frame[1][j + 1 - anti_codon_pos[codon_pos[j]]] != 1 || orfs[i].match_frame[1][j + 1 - anti_codon_pos[codon_pos[j]]] != 1)
				{
					continue;
				}

				int sum_diff = (cod0[0] != cod1[0]) + (cod0[1] != cod1[1]) + (cod0[2] != cod1[2]);
				int anti_sum_diff = (anti_cod0[0] != anti_cod1[0]) + (anti_cod0[1] != anti_cod1[1]) + (anti_cod0[2] != anti_cod1[2]);

				if (sum_diff > 1 || anti_sum_diff>1)
				{
					continue;
				}
				if (!aa_map.count(cod0[0] * 4 * 4 + cod0[1] * 4 + cod0[2]))
				{
					cout << "\nerror: " << cod0[0] << " " << cod0[1] << " " << cod0[2];
					getchar();
				}
				char aa0 = aa_map.at(cod0[0] * 4 * 4 + cod0[1] * 4 + cod0[2]);
				char anti_aa0 = aa_map.at(rc_nuc.at(anti_cod0[2]) * 4 * 4 + rc_nuc.at(anti_cod0[1]) * 4 + rc_nuc.at(anti_cod0[0]));
				for (int k = 0; k < 4; k++)
				{
					if (cod0[h] != k)
					{
						vector<int> pos_cod = cod0;
						pos_cod[h] = k;
						char pos_aa = aa_map.at(pos_cod[0] * 4 * 4 + pos_cod[1] * 4 + pos_cod[2]);
							
						vector<int> pos_anti_cod = anti_cod0;
						pos_anti_cod[anti_codon_pos[h]] = k;
						char pos_anti_aa = aa_map.at(rc_nuc.at(pos_anti_cod[2]) * 4 * 4 + rc_nuc.at(pos_anti_cod[1]) * 4 + rc_nuc.at(pos_anti_cod[0]));
						if ((anti_synonymous && pos_anti_aa != anti_aa0)||(!anti_synonymous && pos_anti_aa==anti_aa0))
						{
							continue;
						}
						if (pos_aa == 'X' || pos_aa == 'M'||pos_anti_aa=='X')
						{
							continue;
						}
						//cout << "\n" << pos_aa << " " << aa0 << " " << pos_anti_aa << " " << anti_aa0 << " " << k << " " << cod0[0]<<cod0[1]<<cod0[2]<<" "<<anti_cod0[0]<<anti_cod0[1]<<anti_cod0[2]<<" "<< cod1[h];
						//getchar();
						//int id = orfs[i].focal_align_seq[1][j - 1] * 4 * 4 * 4 + orfs[i].focal_align_seq[1][j] * 4 * 4 + orfs[i].focal_align_seq[1][j + 1] * 4 + k;
						dnds2.push_back(DNDS_reporter2());
						dnds2.back().context = rev_nuc_map.at(orfs[i].focal_align_seq[1][j - 1]) + rev_nuc_map.at(orfs[i].focal_align_seq[1][j]) + rev_nuc_map.at(orfs[i].focal_align_seq[1][j + 1]);
						dnds2.back().base = rev_nuc_map.at(orfs[i].focal_align_seq[1][j]);
						dnds2.back().sub = rev_nuc_map.at(k);
						dnds2.back().orf_id = orfs[i].orf_id;
						if (pos_aa == aa0)
						{
							//dnds[i].syn_exp += sub_rates[cod0[h]][k];
							//dnds[i].syn_base[cod0[h] * 4 + k]++;
							dnds2.back().syn = 1;
							/////
							
							//cout << "\nsyn: " << dnds2.back().context << " " << dnds2.back().base << " " << dnds2.back().sub << " " << dnds2.back().orf_id<<" h"<<h<<" k"<<k;
							//cout << " " << cod0[0] << cod0[1] << cod0[2] <<" "<<cod1[0]<<cod1[1]<<cod1[2]<<" "<<anti_cod0[0]<<anti_cod0[1]<<anti_cod0[2]<<" "<< anti_cod1[0] << anti_cod1[1] << anti_cod1[2];
							//getchar();

							////
							//dnds2[id].syn_base++;
							if (cod1[h] == k)
							{
								//dnds[i].syn_obs++;
								//dnds[i].syn_subs[cod0[h] * 4 + k]++;
								dnds2.back().actual = 1;
								//dnds2[id].syn_sub++;
							}
						}
						else
						{
							dnds2.back().syn = 0;
							//dnds[i].nonsyn_exp += sub_rates[cod0[h]][k];
							//dnds2[id].nonsyn_base++;
							//dnds[i].nonsyn_base[cod0[h] * 4 + k]++;
							if (cod1[h] == k)
							{
								//dnds[i].nonsyn_obs++;
								//dnds[i].nonsyn_subs[cod0[h] * 4 + k]++;
								dnds2.back().actual = 1;
								//dnds2[id].nonsyn_sub++;
							}
						}
					}
				}

			}
		}
	}
}


void collective_dnds()
{
	vector<string> spp = { "Scer","Spar","Smik","Sjur","Skud","Sarb","Suva","Seub" };
	vector<Orf> orfs;
	read_orfs_long(orfs, "orfs_scer_blast_matched_validated", spp.size());
	cout << "\nread orfs: " << orfs.size();

	//vector<DNDS_reporter> dnds_reporter(orfs.size());
	vector<DNDS_reporter2> dnds_reporter2;

	do_dnds(/*dnds_reporter, */dnds_reporter2, orfs);
	//print_dnds(dnds_reporter, "spar_real_dnds");
	print_dnds2(dnds_reporter2, "spar_dnds_all_");
}

void anti_dnds()
{
	vector<string> spp = { "Scer","Spar","Smik","Sjur","Skud","Sarb","Suva","Seub" };
	vector<Orf> orfs;
	read_orfs_long(orfs, "orfs_scer_blast_matched_validated", spp.size());
	cout << "\nread orfs: " << orfs.size();
	for (int i = 1; i <= 3; i++)
	{
		//vector<DNDS_reporter> dnds_reporter(orfs.size());
		vector<DNDS_reporter2> dnds_reporter2;
		do_anti_dnds(/*dnds_reporter, */dnds_reporter2, orfs, i, true);
		//print_dnds(dnds_reporter, "spar_anti_dnds_antisyn"+to_string(i));
		print_dnds2(dnds_reporter2, "spar_antisyn_all_" + to_string(i));
	}
	for (int i = 1; i <= 3; i++)
	{
		//vector<DNDS_reporter> dnds_reporter(orfs.size());
		vector<DNDS_reporter2> dnds_reporter2;
		do_anti_dnds(/*dnds_reporter, */dnds_reporter2, orfs, i, false);
		//print_dnds(dnds_reporter, "spar_anti_dnds_antinonsyn" + to_string(i));
		print_dnds2(dnds_reporter2, "spar_antinonsyn_all_" + to_string(i));
	}
}

void get_nuc_diverse(vector<Orf> &protogenes, const vector<Vcf> &vcfs, const vector<vector<int>> vcf_map)
{
	for (int i = 0; i < protogenes.size(); i++)
	{
		vector<double> nuc_diverse;
		for (int j = protogenes[i].start_pos; j <= protogenes[i].end_pos; j++)
		{
			if (j >= vcf_map.at(protogenes[i].contig).size())
			{
				cout << "\n" << i << " " << j << " " << protogenes[i].contig << " " << protogenes[i].id;
				getchar();
			}
			int vcf_id = vcf_map.at(protogenes[i].contig).at(j);
			if (vcf_id >= 0 && vcfs[vcf_id].max_length == 1)
			{
				//nuc_diverse.push_back(0);
				vector<int> variant_counts(vcfs[vcf_id].variants.size() + 1);
				vector<double> variant_freqs(vcfs[vcf_id].variants.size() + 1);
				for (int v = 0; v < vcfs.at(vcf_id).calls.size(); v++)
				{
					for (int p = 0; p < 2; p++)
					{
						if (vcfs[vcf_id].calls[v][p] > -1)
						{
							variant_counts[vcfs[vcf_id].calls[v][p]]++;
						}
					}
				}
				int sum = 0;
				for (int n = 0; n < variant_counts.size(); n++)
				{
					sum += variant_counts[n];
				}
				for (int n = 0; n < variant_counts.size(); n++)
				{
					variant_freqs[n] = (double)variant_counts[n] / (double)sum;
				}
				double sumsq = 0;
				for (int n = 0; n < variant_freqs.size(); n++)
				{
					sumsq += variant_freqs[n] * variant_freqs[n];
				}
				//cout << "\n" << i << " " << j << " " << variant_freqs[0] << " " << variant_freqs[1] << " " << sumsq;
				//getchar();
				nuc_diverse.push_back(1 - sumsq);
			}
			else if (vcf_id == -1)
			{
				nuc_diverse.push_back(0);
			}
		}
		double pi = 0;
		for (int n = 0; n < nuc_diverse.size(); n++)
		{
			pi += nuc_diverse[n] / (double)nuc_diverse.size();
		}
		protogenes[i].nuc_diverse = pi;
	}
}

void print_nuc_diverse(vector<Orf> &orfs)
{
	ofstream file("nuc_diverse");
	file<<"id pi";
	for(int i=0;i<orfs.size();i++)
	{
		file<<"\n"<<i<<" "<<orfs[i].nuc_diverse;
	}
}

void nuc_diverse()
{
	vector<Contig> contigs;
	read_fasta(contigs, "S288C_reference_sequence_R64-2-1_20150113.fsa");

	vector<string> spp = { "Scer","Spar","Smik","Sjur","Skud","Sarb","Suva","Seub" };
	vector<Orf> orfs;
	read_orfs_long(orfs, "orfs_scer_blast_matched_validated", spp.size());
	cout<<"\norfs read:" <<orfs.size();
	vector<Vcf> vcfs;
	vector<string> isolates;

	read_vcfs(vcfs, isolates, "1011Matrix.gvcf");
	cout << "\nvcfs read: " << vcfs.size();
	vector<vector<int>> vcf_map;
	map_vcfs_to_genome(vcf_map, vcfs, contigs);
	cout<<"map vcfs";
	
	get_nuc_diverse(orfs, vcfs, vcf_map);
	cout<<"get nuc diverse";
	print_nuc_diverse(orfs);	
}

int main(int argc, char* argv[])
{
	vector<string> args(argv + 1, argv + argc);
	if (args.size() == 0)
	{
		return 0;
	}
	else if (args[0] == "-GetAllORFs")//produce files with information on ORFs in the s cerevisiae genome
	{
		get_all_orfs();
	}
	else if (args[0] == "-MapRiboseqReads")//given a list of fastq files, processes and maps all riboseq reads to genome.
	{
		int start_map = -1;
		int end_map = -1;
		if (args.size()==3)
		{
			start_map = stoi(args[1]);
			end_map = stoi(args[2]);
		}
		map_riboseq_reads(start_map,end_map);
	}
	else if (args[0] == "-CombineRiboseqReads")//assembles and combines mapped riboseq reads from a collection of studies. Must be run after MapRiboseqReads
	{
		int start_map = -1;
		int end_map = -1;
		if (args.size() == 3)
		{
			start_map = stoi(args[1]);
			end_map = stoi(args[2]);
		}
		combine_riboseq_reads(start_map,end_map);
	}
	else if (args[0] == "-IdentifyTranslatedORFs")//from mapped riboseq reads, generates info on riboseq reads associated with each ORF. 
	{//can add an additional flag, described in function: multi_studies chx_studies_sampled ypd_studies ypd_endpoint ypd_no_chx study_accumulation
		string mod = "";
		if (args.size() > 1)
		{
			mod = args[1];
		}
		identify_translated_orfs(mod);
	}
	else if (args[0] == "-GetSubsMatrix")
	{
		get_sub_matrix();
	}
	else if (args[0] == "-SyntenyAlign") // need file: muscle3.8.31_i86linux64
	{
		int species = stoi(args[1]);
		synteny_align(species, false, false);
	}
	else if (args[0] == "-SyntenyMatch")
	{
		synteny_match();
	}
	else if (args[0] == "-SyntenyCheck")
	{
		synteny_check();
	}
	else if (args[0] == "-SensuStrictoBlast")
	{
		sensu_stricto_blast();
	}
	else if (args[0] == "-BlastAlign")
	{
		int species = stoi(args[1]);
		blast_align(species);
	}
	else if (args[0] == "-BlastMatch")
	{
		blast_match();
	}
	else if (args[0] == "-ResolveAlignments")
	{
		resolve_alignments();
	}
	else if (args[0] == "-ValidateAlignments")
	{
		validate_alignments();
	}
	else if (args[0] == "-MultAlignOrthologs")
	{
		mult_align_orthologs();
	}
	else if (args[0] == "-PhyloBlast") // -TBLASTN_SCRAMBLED -ORF_BLAST
	{
		phylo_blast(args[1]);
	}
	else if (args[0] == "-PhyloORFs")
	{
		phylo_find_orfs();
	}
	else if (args[0] == "-PhyloAnalyze") //-TBLASTN_SCRAMBLED
	{
		phylo_analyze(args[1]);
	}
	else if (args[0] == "-SelfBlast")
	{
		self_blast();
	}
	else if(args[0] == "-CodingScores")
	{
		coding_scores();
	}
	else if(args[0] == "-MakeCleanAlignments")
	{
		make_clean_alignment_files();
	}
	else if(args[0] == "-MakeJoinedCleanAlignments")
	{
		make_joined_clean_alignment();
	}
	else if(args[0] == "-IndividualDNDS")
	{
		individual_dnds();
	}
	else if(args[0] == "-StrainSetup")
	{
		strain_setup();
	}
	else if(args[0] == "-StrainAnalyze")
	{
		strain_analyze();
	}
	else if(args[0] == "-NucDiverse")
	{
		nuc_diverse();
	}
	else if(args[0] == "-CollectiveDNDS")
	{
		collective_dnds();
	}
	else if(args[0] == "-AntiDNDS")
	{
		anti_dnds();
	}
}