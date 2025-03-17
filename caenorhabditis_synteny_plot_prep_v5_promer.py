# caenorhabditis_synteny_plot_prep.py

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys, os, getopt, datetime

class Gene:
    def __init__(self, gene_id, contig, start, end, orientation):
        self.gene_id = gene_id
        self.contig = contig
        self.start = start
        self.end = end
        self.orientation = orientation
    '''def get_contig(self):
        return self.contig
    def get_start(self):
        return(self.start)
    def get_end(self):
        return self.end
    def get_id(self):
        return self.gene_id
    '''

# round up to the closest 1000 nucleotide coordinate
def round_up_1000(coordinate):
    return coordinate + 1000 - coordinate%1000

# round down to the closest 1000 nucleotide coordinate plus 1
def round_down_1000(coordinate):
    return coordinate + 1 - coordinate%1000

# C elegans and C briggsae alternate transcripts have "a", "b", "c" suffixes
def a_isoform_only(input_protein_fasta, output_protein_fasta):
    protein_list = []
    with open(input_protein_fasta, "r") as in_fh, open(output_protein_fasta, "w") as out_fh:
        for protein_sequence in SeqIO.parse(in_fh, "fasta"):
            if protein_sequence.id[-1].isalpha() == True:
                if protein_sequence.id[-1] == "a":
                    SeqIO.write([protein_sequence], out_fh, "fasta")
                    protein_list.append(protein_sequence.id)
                else:
                    pass
            else:
                SeqIO.write([protein_sequence], out_fh, "fasta")
                protein_list.append(protein_sequence.id)
    return protein_list

# others have "t1", "t2", "t3" transcripts
def t1_isoforms_only(input_protein_fasta, output_protein_fasta):
    protein_list = []
    with open(input_protein_fasta, "r") as in_fh, open(output_protein_fasta, "w") as out_fh:
        for protein_sequence in SeqIO.parse(in_fh, "fasta"):
            id_length = len(protein_sequence.id)
            if protein_sequence.id[id_length-2:id_length] == "t1":
                SeqIO.write([protein_sequence], out_fh, "fasta")
                protein_list.append(protein_sequence.id)
            else:
                pass
    return protein_list

# get a list of non-redundant protein names from the protein_fasta
def get_non_redundant_protein_names(protein_fasta):
    protein_list = []
    for protein_sequence in SeqIO.parse(protein_fasta, "fasta"):
        protein_list.append(protein_sequence.id)
    return protein_list

def get_gff_attribute(gff_line, attribute):
    attributes = gff_line.split("\t")[8].rstrip()
    
    attributes_dict = {}
    # I	WormBase	mRNA	11495	16837	.	+	.	ID=Transcript:Y74C9A.2a.1;Parent=Gene:WBGene00022276;Name=Y74C9A.2a.1;wormpep=CE24660;locus=nlp-40;uniprot_id=Q9N4D8
    for item in attributes.split(";"):
        attributes_dict[item.split("=")[0]] = item.split("=")[1]
    
    if attribute in attributes_dict:
        return attributes_dict[attribute]
    else:
        return "NA"

def get_gff_start(gff_line):
    return int(gff_line.split("\t")[3])

def get_gff_end(gff_line):
    return int(gff_line.split("\t")[4])
    
def get_gff_orientation(gff_line):
    return gff_line.split("\t")[6]

def get_gff_type(gff_line):
    return gff_line.split("\t")[2]
    
def get_gff_contig(gff_line):
    return gff_line.split("\t")[0]


# read_gff
# returns:
#   - a gene_dict[gene_id] that can return the contig on which the gene resides and the gene number index
#   - a contig_dict of lists pointing to Gene objects
# suffix: some gffs use ".1" after the protein name to indicate mRNA
#   eg. protein = Y74C9A.2a, mRNA_Name = Y74C9A.2a.1
def read_gff(gff_file, suffix, species_abbrev, protein_list):
    contig_dict = {}
    gene_dict = {}
    with open(gff_file, "r") as gff_fh:
        line = gff_fh.readline()
        
        while line:
            if line[0] != "#":
                if get_gff_type(line) == "mRNA":
                    gene_id = get_gff_attribute(line, "Name")
                        
                    if suffix == ".1":
                        gene_id = gene_id[0:len(gene_id)-2]
                        
                        
                    if gene_id in protein_list:
                        contig = species_abbrev + "." + get_gff_contig(line)
                        start = get_gff_start(line)
                        end = get_gff_end(line)
                        orientation = get_gff_orientation(line)
                            
                        new_gene = Gene(gene_id, contig, start, end, orientation)
                            
                        if contig not in contig_dict:
                            contig_dict[contig] = [new_gene]
                        else:
                            contig_dict[contig].append(new_gene)
                    else:
                        pass
                else:
                    pass
            else:
                pass
            line = gff_fh.readline()
    # sort the lists of genes on each contig by starting position
    for contig in contig_dict:
        contig_dict[contig].sort(key=lambda x: x.start)
        
        # now populate the gene_dict
        # key = gene_id
        # value = contig, index
        # such that
        # gene_dict[gene_id] = [contig, index]
        for index, gene_object in enumerate(contig_dict[contig]):
            # index is the index at which gene_object is found
            gene_dict[gene_object.gene_id] = [contig, index]
    
    return contig_dict, gene_dict
    
# get write_contig_block_to_file(fasta_file, contig_name, start, end)
def write_contig_block_to_file(fasta_file, species_abbrev, contig_name, start, end, working_directory):
    temp_contig_name = contig_name.replace(species_abbrev + ".", "")
    
    new_seq = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))[temp_contig_name].seq
    
    start = max([1, start])
    end = min([len(new_seq), end])
    
    new_seq = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))[temp_contig_name].seq[start-1:end]
    new_id = contig_name + "." + str(start) + "-" + str(end)
    new_fasta_file = working_directory + "/" + new_id + ".DNA.fasta"
    
    new_seqrecord = SeqRecord(new_seq, id=new_id, description="")
    
    SeqIO.write([new_seqrecord], new_fasta_file, "fasta")
    
    # return the newly created path:file 
    return new_fasta_file
    
    
            
# config file should have columns like this:
# species	abbrev	raw_proteins_file	curated_proteins_file	gff_file    genome_file
def read_config_file(config_file, raw_proteins_dict, curated_proteins_dict, gff_dict, genome_dict):
    species_abbrev_list = []
    with open(config_file, "r") as config_fh:
        line = config_fh.readline()
        
        while line:
            if line[0] == "#":
                pass
            else:
                species_abbrev = line.split("\t")[1]
                species_abbrev_list.append(species_abbrev)
                raw_proteins_dict[species_abbrev] = line.split("\t")[2]
                curated_proteins_dict[species_abbrev] = line.split("\t")[3]
                gff_dict[species_abbrev] = line.split("\t")[4]
                genome_dict[species_abbrev] = line.split("\t")[5].rstrip()
            line = config_fh.readline()
    return species_abbrev_list
            


def main(argv):
    arg_species = "Celeg"
    arg_gene = "C38D4.4"
    arg_e = "0.001"
    arg_block_size = 5
    arg_l = "6"
    arg_c = "20"
    working_directory = ""
    config_file = "caenorhabditis_synteny_config.tsv"
    arg_help = "{0} [options] -s <species> -g <gene> -e <mmseqs_evalue> -b <synteny_block_size> -l <promer_min_match_length> -c <promer_min_cluster> -p <promer_only>\n".format(argv[0])
    
    # get the options
    opts, args = getopt.getopt(argv[1:], "hs:g:e:b:l:c:p:C:", ["help", "species=", "gene=", "mmseqs_evalue=", "synteny_block_size=", "promer_min_match_length=", "promer_min_cluster=", "promer_only_directory=", "config="])
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(0)
        elif opt in ("-s", "--species"): 
            arg_species = arg
        elif opt in ("-g", "--gene"):
            arg_gene = arg
        elif opt in ("-e", "--mmseqs_evalue"):
            arg_e = arg
        elif opt in ("-b", "--synteny_block_size"):
            arg_block_size = int(arg)
        elif opt in ("-l", "--promer_min_match_length"):
            arg_l = arg
        elif opt in ("-c", "--promer_min_cluster="):
            arg_c = arg
        elif opt in ("-p", "--promer_only_directory="):
            working_directory = arg
        elif opt in ("-C", "--config="):
            config_file = arg
    
    # declare block_list
    # a list of FASTA files containing the putative syntenic blocks
    block_list = []
    
    if working_directory == "":
        #######################################################################################################################################################    
        #######################################################################################################################################################
        # PREPARE THE PROTEIN FASTA FILES AND READ FROM THE GFF3 FILES
        #######################################################################################################################################################
        #######################################################################################################################################################
        
        # make a working directory
        working_directory = str(datetime.datetime.now()).replace(" ", "_").replace(":", "") + "." + arg_species + "." + arg_gene
        os.system("mkdir " + working_directory)
        
        raw_proteins_dict = {}           
        curated_proteins_dict = {}
        gff_dict = {}
        genome_dict = {}
        
        species_abbrev_list = read_config_file(config_file, raw_proteins_dict, curated_proteins_dict, gff_dict, genome_dict)
        
        # capture the run command in a log file "log.txt"
        log_file = f"{working_directory}/log.txt"
        with open(log_file, "w") as log_fh:
            for item in argv:
                log_fh.write(item + " ")
            log_fh.write("\n\narguments:\n")
            log_fh.write(f"-s = {arg_species}, query species\n")
            log_fh.write(f"-g = {arg_gene}, query gene\n")
            log_fh.write(f"-b = {arg_block_size}, number of genes before and after the query\n")
            log_fh.write(f"-e = {arg_e}, mmseqs e-value\n")
            log_fh.write(f"-l = {arg_l}, promer -l argument: minimum match length\n")
            log_fh.write(f"-c = {arg_c}, promer -c argument: minimum cluster size\n")
            log_fh.write(f"-C = {config_file}, configuration file listing protein, gff3 and genome files\n")
            log_fh.write(f"available species: ")
            for index, species_abbrev in enumerate(species_abbrev_list):
                log_fh.write(species_abbrev)
                if index < len(species_abbrev_list) - 1:
                    log_fh.write(", ")
                else:
                    log_fh.write("\n\n")
                
        
        # protein_list
        # protein_list["Celeg"] = ["gene1", "gene2", "gene3", ...]
        protein_list = {}
        
        # contig_dict
        # contig_dict[species][contig][index] = Gene_object
        contig_dict = {}
        gene_dict = {}
        
        for species_abbrev in species_abbrev_list:
            a_isoform_list = ["Celeg",
                              "Cbrig",
                              "Cbren"]
            if os.path.exists(curated_proteins_dict[species_abbrev]) == False:
                if species_abbrev in a_isoform_list:
                    protein_list[species_abbrev] = a_isoform_only(raw_proteins_dict[species_abbrev], curated_proteins_dict[species_abbrev])
                else:
                    protein_list[species_abbrev] = t1_isoforms_only(raw_proteins_dict[species_abbrev], curated_proteins_dict[species_abbrev])
            else:
                protein_list[species_abbrev] = get_non_redundant_protein_names(curated_proteins_dict[species_abbrev])
            
            # capture the ".1" mRNA naming scheme
            the_point_one_mRNA_list = ["Celeg",
                                       "Cbrig",
                                       "Cbren"]
            
            suffix = ""
            if species_abbrev in the_point_one_mRNA_list:
                suffix = ".1"
            else:
                pass
            
            contig_dict[species_abbrev], gene_dict[species_abbrev] = read_gff(gff_dict[species_abbrev], suffix, species_abbrev, protein_list[species_abbrev])
            
            '''
            with open(working_directory + "/" + species_abbrev + ".contig_dict.txt", "w") as contig_fh, open(working_directory + "/" + species_abbrev + ".gene_dict.txt", "w") as gene_fh:
                for contig in contig_dict[species_abbrev]:
                    for gene_item in contig_dict[species_abbrev][contig]:
                        line_to_write = contig + "\t" + gene_item.gene_id + "\t" + str(gene_item.start) + "\t" + str(gene_item.end) + "\t" + gene_item.orientation + "\n"
                        contig_fh.write(line_to_write)
                # recall that gene_dict[gene_id] = [contig, index]
                for gene_name in gene_dict[species_abbrev]:
                    contig = gene_dict[species_abbrev][gene_name][0]
                    index = gene_dict[species_abbrev][gene_name][1]
                    line_to_write = gene_name + "\t" + contig + "\t" + str(index) + "\n"
                    gene_fh.write(line_to_write)
            '''       
        #######################################################################################################################################################    
        #######################################################################################################################################################
        # IDENTIFY PUTATIVE SYNTENY BLOCKS
        #######################################################################################################################################################
        #######################################################################################################################################################
        
        gggenomes_gene_file = working_directory + "/gggenomes_gene_file.txt"
        gggenomes_seq_file = working_directory + "/gggenomes_seq_file.txt"
        
        
        
        with open(gggenomes_gene_file, "w") as gggenomes_gene_fh, open(gggenomes_seq_file, "w") as gggenomes_seq_fh:
            gggenomes_gene_fh.write("seq_id\tstart\tend\tstrand\tfeat_id\twidth\tname\tgeom_id\ttype\n")
            gggenomes_seq_fh.write("seq_id\tlength\n")
            
            # grab the query gene as well and write to file
            [query_contig, query_index] = gene_dict[arg_species][arg_gene]
            query_gene_fasta = working_directory + "/" + query_contig + "." + arg_gene + ".fasta"
            query_seqrecord = SeqIO.to_dict(SeqIO.parse(curated_proteins_dict[arg_species], "fasta"))[arg_gene]
            SeqIO.write([query_seqrecord], query_gene_fasta, "fasta")
            
            # write down all putative orthologues of the query gene
            putative_orthologues_fasta = f"{working_directory}/{arg_species}.{arg_gene}.putativeOrthologues.fasta"
            
            with open(putative_orthologues_fasta, "w") as ortho_fh:
                SeqIO.write([query_seqrecord], ortho_fh, "fasta")
                
                # now go through all the species and grab the potential syntenic block
                for species_abbrev in species_abbrev_list:
                    if species_abbrev == arg_species:
                        # use the query species_abbrev, and gene
                        # and grab the DNA that contains the gene +/- arg_block_size genes
                        # grab the block of genes (eg. 5 genes before and after the query gene)
                        block_start_index = max([query_index-arg_block_size, 0])
                        block_end_index = min([query_index+arg_block_size, len(contig_dict[arg_species][query_contig])-1])
                        
                        # grab the absolute coordinates, round to the nearest 1000s
                        block_start = round_down_1000(contig_dict[arg_species][query_contig][block_start_index].start)
                        block_end = round_up_1000(contig_dict[arg_species][query_contig][block_end_index].end)
                        
                        # write down the gene coordinates
                        for index in range(block_start_index, block_end_index+1):
                            seq_id = query_contig + "." + str(block_start) + "-" + str(block_end)
                            start = str(contig_dict[arg_species][query_contig][index].start - block_start + 1)
                            end = str(contig_dict[arg_species][query_contig][index].end - block_start + 1)
                            strand = contig_dict[arg_species][query_contig][index].orientation
                            feat_id = contig_dict[arg_species][query_contig][index].gene_id
                            width = contig_dict[arg_species][query_contig][index].end - contig_dict[arg_species][query_contig][index].start + 1
                            name = feat_id
                            geom_id = feat_id
                            gggenomes_gene_fh.write(seq_id+"\t"+str(start)+"\t"+str(end)+"\t"+strand+"\t"+feat_id+"\t"+str(width)+"\t"+name+"\t"+geom_id+"\tCDS\n") 
                        
                        query_DNA_block = write_contig_block_to_file(genome_dict[arg_species], arg_species, query_contig, block_start, block_end, working_directory)
                        gggenomes_seq_fh.write(query_contig + "." + str(block_start) + "-" + str(block_end)  + "\t" + str(block_end - block_start + 1) + "\n")
                        
                        # record down this file so we may do promer with it
                        block_list.append(query_DNA_block)

                    
                    else:
                        # do an mmseqs search for the likely ortholog
                        temp_search_output = "temp.mmseqs.output.txt"
                        temp_folder = "temp"
                        temp_error_messages = temp_search_output + ".err"
                        line_to_write = "mmseqs easy-search -e " + arg_e + " " + query_gene_fasta + " " + curated_proteins_dict[species_abbrev] + " " + temp_search_output + " " + temp_folder + " > " + temp_error_messages + " 2>&1"
                        os.system(line_to_write)
                                
                        # fish out the mmseqs results
                        hit_dict = SeqIO.to_dict(SeqIO.parse(curated_proteins_dict[species_abbrev], "fasta"))
                        with open(log_file, "a") as log_fh:
                            log_fh.write("\nreading " + curated_proteins_dict[species_abbrev] + "\n")
                        
                        with open(temp_search_output, "r") as temp_search_output_fh:
                            line = temp_search_output_fh.readline()
                                    
                            while line:
                                hit_gene_id = line.split("\t")[1]
                                with open(log_file, "a") as log_fh:
                                    log_fh.write(hit_gene_id + "\n")
                                
                                SeqIO.write([hit_dict[hit_gene_id]], ortho_fh, "fasta")
                                
                                # grab the block of genes including the hit (eg. 5 genes before and after the hit gene)
                                [hit_contig, hit_index] = gene_dict[species_abbrev][hit_gene_id]
                                block_start_index = max([hit_index-arg_block_size, 0])
                                block_end_index = min([hit_index+arg_block_size, len(contig_dict[species_abbrev][hit_contig])-1])
                                
                                block_start = round_down_1000(contig_dict[species_abbrev][hit_contig][block_start_index].start)
                                block_end = round_up_1000(contig_dict[species_abbrev][hit_contig][block_end_index].end)
                                
                                # write down the gene coordinates
                                for index in range(block_start_index, block_end_index+1):
                                    seq_id = hit_contig + "." + str(block_start) + "-" + str(block_end)
                                    start = str(contig_dict[species_abbrev][hit_contig][index].start - block_start + 1)
                                    end = str(contig_dict[species_abbrev][hit_contig][index].end - block_start + 1)
                                    strand = contig_dict[species_abbrev][hit_contig][index].orientation
                                    feat_id = contig_dict[species_abbrev][hit_contig][index].gene_id
                                    width = contig_dict[species_abbrev][hit_contig][index].end - contig_dict[species_abbrev][hit_contig][index].start + 1
                                    name = feat_id
                                    geom_id = feat_id
                                    gggenomes_gene_fh.write(seq_id+"\t"+str(start)+"\t"+str(end)+"\t"+strand+"\t"+feat_id+"\t"+str(width)+"\t"+name+"\t"+geom_id+"\tCDS\n") 
                                
                                hit_DNA_block = write_contig_block_to_file(genome_dict[species_abbrev], species_abbrev, hit_contig, block_start, block_end, working_directory)
                                
                                if hit_DNA_block not in block_list:
                                    block_list.append(hit_DNA_block)
                                    gggenomes_seq_fh.write(hit_contig + "." + str(block_start) + "-" + str(block_end) + "\t" + str(block_end - block_start + 1) + "\n")
                                else:
                                    pass
                                
                                line = temp_search_output_fh.readline()
    
    else:
        for entry in os.scandir(working_directory):
            if ".DNA.fasta" in entry.path:
                block_list.append(entry.path)
            else:
                pass
    #######################################################################################################################################################    
    #######################################################################################################################################################
    # USE PROMER TO IDENTIFY SYNTENIC STRETCHES OF DNA
    #######################################################################################################################################################
    #######################################################################################################################################################    
    # do promer, two species at a time
    gggenomes_link_file = f"{working_directory}/gggenomes_link_file.txt"
    
    with open(gggenomes_link_file, "w") as link_fh:
        link_fh.write("seq_id\tstart\tend\tseq_id2\tstart2\tend2\n")
        for block_A in block_list:
            for block_B in block_list:
                if block_A != block_B:
                    os.system(f"promer -c {arg_c} -l {arg_l} -p temp {block_A} {block_B}")
                    os.system("show-coords -T temp.delta > temp.coords")
                        
                    with open("temp.coords", "r") as coords_fh:
                        for i in range(4):
                            coords_fh.readline()
                        line = coords_fh.readline()
                        while line:
                            line_items = line.rstrip().split("\t")
                            seq_id = line_items[11]
                            start = line_items[0]
                            end = line_items[1]
                            seq_id2 = line_items[12]
                            start2 = line_items[2]
                            end2 = line_items[3]
                            link_fh.write(seq_id + "\t" + start + "\t" + end + "\t" + seq_id2 + "\t" + start2 + "\t" + end2 + "\n")
                            line = coords_fh.readline()
                        
                    os.system(f"promer -c {arg_c} -l {arg_l} -p temp {block_B} {block_A}")
                    os.system("show-coords -T temp.delta > temp.coords")
                        
                    with open("temp.coords", "r") as coords_fh:
                        for i in range(4):
                            coords_fh.readline()
                        line = coords_fh.readline()
                        while line:
                            line_items = line.rstrip().split("\t")
                            seq_id = line_items[11]
                            start = line_items[0]
                            end = line_items[1]
                            seq_id2 = line_items[12]
                            start2 = line_items[2]
                            end2 = line_items[3]
                            link_fh.write(seq_id + "\t" + start + "\t" + end + "\t" + seq_id2 + "\t" + start2 + "\t" + end2 + "\n")
                            line = coords_fh.readline()
    os.system("rm -r temp*")
    
    #######################################################################################################################################################    
    #######################################################################################################################################################
    # EVLUATE WHICH SYNTENY BLOCK IS THE BEST TO GRAPH
    #######################################################################################################################################################
    #######################################################################################################################################################
    
    # link_arrays[blockA][blockB] = [0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0]
    # link_arrays[blockA]["length"] = 29398479
    link_arrays = {}
    
    for block_A in block_list:
        block_A_name = SeqIO.read(block_A, "fasta").id
        block_A_length = len(SeqIO.read(block_A, "fasta"))
        link_arrays[block_A_name] = {"length": block_A_length}
        for block_B in block_list:
            if block_B != block_A:
                block_B_name = SeqIO.read(block_B, "fasta").id
                #print(block_A_name)
                #print(block_B_name, "\n")
                link_arrays[block_A_name][block_B_name] = [0]*(block_A_length + 1)
            else:
                pass
                
    with open(gggenomes_link_file, "r") as link_fh:
        link_fh.readline() # header line
        
        line = link_fh.readline()
        
        while line:
            block_A_name = line.split("\t")[0]
            block_A_start = int(line.split("\t")[1])
            block_A_end = int(line.split("\t")[2])
            block_B_name = line.split("\t")[3]
            block_B_start = int(line.split("\t")[4])
            block_B_end = int(line.rstrip().split("\t")[5])
            
            for i in range(block_A_start, block_A_end + 1):
                #print(block_A_name, " ", block_A_start, block_B_name, " ", block_B_start)
                link_arrays[block_A_name][block_B_name][i] = 1
            
            for i in range(block_B_start, block_B_end + 1):
                link_arrays[block_B_name][block_A_name][i] = 1
            
            line = link_fh.readline()
        
    total_links_file = f"{working_directory}/total_links_file.tsv"
    
    with open(total_links_file, "w") as total_links_fh:
        for index, block_B in enumerate(block_list):
            block_B_name = SeqIO.read(block_B, "fasta").id
            if index == 0:
                total_links_fh.write("\t" + block_B_name + "\t")
            elif index == len(block_list) - 1:
                total_links_fh.write(block_B_name + "\n")
            else:
                total_links_fh.write(block_B_name + "\t")
                
        for block_A in block_list:
            block_A_name = SeqIO.read(block_A, "fasta").id
            total_links_fh.write(block_A_name)
            
            for index, block_B in enumerate(block_list):
                block_B_name = SeqIO.read(block_B, "fasta").id
                
                if block_A_name == block_B_name:
                    total_links_fh.write("\t0")
                else:
                    links_fraction = sum(link_arrays[block_A_name][block_B_name])/link_arrays[block_A_name]["length"]
                    total_links_fh.write("\t" + str(links_fraction))
                
                if index == len(block_list) - 1:
                    total_links_fh.write("\n")
                else:
                    pass
    
    
if __name__ == "__main__":
    main(sys.argv)