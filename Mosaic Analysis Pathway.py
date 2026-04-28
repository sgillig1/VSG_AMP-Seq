import argparse # to take the arguments into the file
from datetime import datetime # time specific tests for optimization
from sys import argv 
import os 
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import csv
import pysam
import time
from Bio import SeqIO

import HTSeq
import sys
import logging


#from Bio.SeqIO.QualityIO import FastqGeneralIterator
from difflib import SequenceMatcher
    

def VSG(reads_R1, reads_R2, reference, final_sam, genome_name, gff_file, target_gene):
    # reads_R2 = fq file name Read 2 (read from 221 primer)
    # reads_R1 = fq file name Read 1 (opposite end read)
    # reference = reference file name (The AM_VSG genome) (Base name for the reference genome)
    # final_sam = output file base name (no extension)
    # genome_name fasta file
    # gff_file = gff file
    # target_gene = target gene id
    
    print("starting " + final_sam + "........")
    reference_VSG = [] #list of all the VSG genes
    
    for record in SeqIO.parse(genome_name, "fasta"):
        reference_VSG.append(record.id)
    #print(reference_VSG)
    
    bowtie_code = "bowtie -p 6 -v 1 -m 1 -X 1500 --quiet --un %s -S %s %s > %s"
    seq_length = 0 # to ensure that the final read doesn't get run multiple times (only if reads same length). This is used when the number of reads is very small

    remove_list = []
    sam_list = [] # list of sam file names to allow merging in the trimmer (this allows for the creation of many files that can be put together)
    for n in range(105):
        sam_list.append(final_sam + str(n) + ".sam")
    new_list = [] #list of the new merged files
    for m in range(105):
        new_list.append(final_sam + "new_" + str(m) + ".sam")
        
    aligned_list = [] #list of the new merged files
    for m in range(105):
        aligned_list.append(final_sam + "aligned_" + str(m) + ".sam")
    
    output_readnames = open("%s_output_readnames.fq" % (final_sam), 'w')
    junk = open("%s_junk.fq" % (final_sam), 'w') # JUNK should contain a list of all the reads from the trimmer that did not end up aligning
    junk_name = final_sam + "_junk.fq"
    
    ##########################################################################
    
    def trimmer(reads): # outputs the sam file name of the trimmed reads, imputs an fq file name of unaligned reads
        i = 1
        junk = open("%s_junk.fq" % (final_sam), 'w') # JUNK should contain a list of all the reads from the trimmer that did not end up aligning
        junk_name = final_sam + "_junk.fq"
        
        final_unaligned = final_sam + "unaligned.fq" # the unaligned file name for all
        while i <= 100: #iterate through the trimming process 100 times
                
            if i == 1: # run first alignment(There is nothing for it to merge with for the aligned SAM) (can you only output the aligned to SAM and not the unaligned??)
                os.system(bowtie_code % (final_unaligned, reference, reads, sam_list[i])) # GET RID OF THIS
                
            elif i == 2: # run alignment and merge the two sam files, delete the first SAM file
                os.system(bowtie_code % (final_unaligned, reference, reads, sam_list[i]))
                os.system("samtools merge -f %s %s %s" % (new_list[i], sam_list[i-1], sam_list[i])) # merge the first two sam files
                os.system("samtools view -F 4 -h %s > %s" % (new_list[i], aligned_list[i])) # only take aligned files
                os.remove(sam_list[i-1]) # delete first SAM
                os.remove(sam_list[i])
            else: # run alignment and merge the two sam files, delete the previous SAM file and the previous merge file
                os.system(bowtie_code % (final_unaligned, reference, reads, sam_list[i]))
                os.system("samtools merge -f %s %s %s" % (new_list[i], aligned_list[i-1], sam_list[i]))
                os.system("samtools view -F 4 -h %s > %s" % (new_list[i], aligned_list[i])) # only take aligned files
                os.remove(new_list[i]) # delete the old new file
                os.remove(aligned_list[i-1]) # delete the old new file
                os.remove(sam_list[i]) # delete the last file
                
            
            outfile_name = final_sam + "trimmed_unaligned.fq" #new file name for trimmed file
            unaligned_reads = open(final_unaligned, 'rU') # open the unaligned reads from bowtie
        
            outfile = open(outfile_name, 'w')
            outfile.truncate(0) #clear the new file
            
            
            for title, seq, qual in FastqGeneralIterator(unaligned_reads): #iterate through the fq entries
                seq_length = len(seq) - 1 # new size of sequence (one bp smaller) this is for the next iteration
                #forwards reads
                if seq_length <= 20:
                    junk.write('@' + title + '\n')
                    junk.write(seq[:seq_length] + '\n')
                    junk.write('+' + '\n')
                    junk.write(qual[:seq_length] + '\n')
                elif i == 100: # Take the junk reads from after trimming
                    if title[-1] != 'b': # if it is a forward read, 
                    
                        junk.write('@' + title + '\n')
                        junk.write(seq[:seq_length] + '\n')
                        junk.write('+' + '\n')
                        junk.write(qual[:seq_length] + '\n')
    
                    elif title[-1] == 'b': # if backwards read 
                        junk.write('@' + title + '\n')
                        junk.write(seq[1:] + '\n')
                        junk.write('+' + '\n')
                        junk.write(qual[1:] + '\n')
                else:
                    if i == 1 or title[-1] != 'b': # if it is a forward read, cut off from the end. (or if in the first iteration)
                    
                        outfile.write('@' + title + '\n')
                        outfile.write(seq[:seq_length] + '\n')
                        outfile.write('+' + '\n')
                        outfile.write(qual[:seq_length] + '\n')
                    
                    #backwards reads 
                    if i == 1: # add the b to the backwards reads for the first iteration
                        outfile.write('@' + title + 'b' + '\n')
                        outfile.write(seq[1:] + '\n')
                        outfile.write('+' + '\n')
                        outfile.write(qual[1:] + '\n')
                    elif title[-1] == 'b': # if backwards read then cut from beginning
                        outfile.write('@' + title + '\n')
                        outfile.write(seq[1:] + '\n')
                        outfile.write('+' + '\n')
                        outfile.write(qual[1:] + '\n')
                
            outfile.close()
            reads = outfile_name # the reads analyzed in next rotation are the trimmed reads
            sam_filename = aligned_list[i] #the current full list of alignments
            
            i = i + 1
        os.remove(new_list[2]) # delete the old new file
        
    
        return sam_filename #returns the file name that contains all the aligned
    
    def locations_trim(reads): #input the sam file from the trimmer
        sam = reads
        bam = sam.split(".")[0] + ".bam"
        sorted_bam = bam.split(".")[0] + ".sorted.bam"
        
        remove_list.append(bam)
        remove_list.append(sorted_bam)
        remove_list.append(sorted_bam + ".bai")
        
        os.system("samtools view -bS -o %s %s" % (bam, sam)) # Turn the sam file into a sorted BAM with an index
        os.system("samtools sort %s -o %s" % (bam, sorted_bam))
        os.system("samtools index %s" % (sorted_bam))
        
        bamfile = pysam.AlignmentFile(sorted_bam, "rb") #uses a bam file with an index
        align_221 = [] # all the reads that align to 221
        new_reads = [] # the opposite reads that corresponds to 221
        positions = [] # the left most base pair of the aligned reads + the length of the read (********** HOW TO CHANGE TO SHOW THAT VSG CAN ALIGN ON EITHER SIDE OF 221? *********)
        read_mosaic = [] # what the other side maps to
        new_positions = [] # the position on 221 that corresponds to read_mosaic
        new_name = [] # the read names
        
        alignment_dict = {}
        reverse_dict = {}
        
        for read in bamfile.fetch(target_gene): # find all the reads that align to 221
            #align_221.append(str(read.query_name))  # list of the aligned reads
            #positions.append(str(read.pos + len(read.query_alignment_sequence))) # pos is the left most position and the len is the legth of the read (this shows the right most base pair in the read)
            
            alignment_dict[str(read.query_name)] = str(read.pos + len(read.query_alignment_sequence))
            
            # if read[-1] == 'b':
            #     positions.append(str(read.pos))
            # else:
            #     positions.append(str(read.pos + len(read.query_alignment_sequence)))
        
        for i in alignment_dict: #retrieve the read that corresponds to the trimmed sequence align to 221
            if i[-1] == 'b':
                reverse_dict[i[:-1]] = [alignment_dict[i]]
            else:
                reverse_dict[i + 'b'] = [alignment_dict[i]]
    
        
        # j = 0
        # for i in align_221: #retrieve the read that corresponds to the trimmed sequence align to 221
        #     if (align_221[j])[-1] == 'b':
        #         new_reads.append((align_221[j])[:-1])
        #     else:
        #         new_reads.append(align_221[j] + 'b')
        #     j = j + 1
        #         
        # j = 0
        
        indexed_file = pysam.IndexedReads(bamfile) # make an indexed read file
        indexed_file.build()
        
        k = -1 #since it doesn't always reach the end of the loop
        for name in reverse_dict:
            k = k + 1
            try:
                iterator = indexed_file.find(name) # check to see if there is a read with the name
            except KeyError:
                pass
            else:    
                for x in iterator: # x is an alignment
                    try:
                        test = str(x.reference_name) #check to see if there is a match/ alignment for the read
                    except Exception as e:
                        print (e)
                        pass
                    else:
                        if str(x.reference_name) != target_gene: # we don't want any 221
                            reverse_dict[name].append(str(x.reference_name)) # What does the other end of the read map to
                            reverse_dict[name].append(str(x.pos))
        return reverse_dict # get into a file and then visualize (histogram like?)
    
    def locations_gap(reads, case2): # input the sam file and the case 2 fq file
        sam = reads
        bam = sam.split(".")[0] + ".bam"
        sorted_bam = bam.split(".")[0] + ".sorted.bam"
        
        remove_list.append(bam)
        remove_list.append(sorted_bam)
        remove_list.append(sorted_bam + ".bai")
        
        os.system("samtools view -bS -o %s %s" % (bam, sam)) # Turn the sam file into a sorted BAM with an index
        os.system("samtools sort %s -o %s" % (bam, sorted_bam))
        os.system("samtools index %s" % (sorted_bam))
        
        bamfile = pysam.AlignmentFile(sorted_bam, "rb") #uses a bam file with an index
        mosaics = {} # all the reads that align to 221
        output = {}
    
        
        for reference in reference_VSG:
            if reference != target_gene:
                for read in bamfile.fetch(reference): # find all the reads that align to 221
                    mosaics[str(read.query_name)] = reference
            
        
        sam = case2
        bam = sam.split(".")[0] + ".bam"
        sorted_bam = bam.split(".")[0] + ".sorted.bam"
        
        remove_list.append(bam)
        remove_list.append(sorted_bam)
        remove_list.append(sorted_bam + ".bai")
        
        os.system("samtools view -bS -o %s %s" % (bam, sam)) # Turn the sam file into a sorted BAM with an index
        os.system("samtools sort %s -o %s" % (bam, sorted_bam))
        os.system("samtools index %s" % (sorted_bam))
        
        bamfile = pysam.AlignmentFile(sorted_bam, "rb") #uses a bam file with an index
        
        k = -1 #since it doesn't always reach the end of the loop
        indexed_file = pysam.IndexedReads(bamfile) # make an indexed read file
        indexed_file.build()
        for query in mosaics:
            original_query = query
            query = query[:-1] + "2"
            #print query
            k = k + 1
            try:
                iterator = indexed_file.find(query) # check to see if there is a read with the name
            except KeyError:
                pass
            else:    
                for x in iterator: # x is an alignment
                    try:
                        test = str(x.reference_name) #check to see if there is a match/ alignment for the read
                    except Exception as e:
                        print (e)
                        pass
                    else:
                        #output[mosaics[query]].append(str(x.pos + len(x.query_alignment_sequence)))
                        #print (str(x.pos + len(x.query_alignment_sequence)))
                        output[original_query] = [mosaics[original_query], "after " + str(x.pos + len(x.query_alignment_sequence))]
        return output
    
    def remove_files(remove):
        print ("Removing files .......")
        remove = list(set(remove))
        for x in remove:
            #print x
            os.remove(x)
    
    def count_reads(filename):
        i = 0
        for title in FastqGeneralIterator(filename):
            i = i+1
        validate_count = i
        return validate_count

    
    ##########################################################################
    
    filename = open(reads_R2, 'rU')
    validate_count = count_reads(filename)
    
    # Align the first read
    unaligned_read1 = final_sam + "unaligned_read1.fq" # unaligned file for first read
    aligned_read1 = final_sam + "aligned_read1.sam" # aligned file name for read 1
    
    os.system(bowtie_code % (unaligned_read1, reference, reads_R2, aligned_read1))
    
    start_time = time.time()
    case1 = trimmer(unaligned_read1) # run trimmer on the unaligned reads
    
    ###########################  SAVE ALL TRIMMED READ FILES AS BAM
    sam_file = case1
    bam_file = final_sam + "Case1_aligned_trim.bam"
    sorted_bam = bam_file.split(".")[0] + ".sorted.bam"
    
    os.system("samtools view -bS -o %s %s" % (bam_file, sam_file))
    os.system("samtools sort %s -o %s" % (bam_file, sorted_bam))
    os.system("samtools index %s" % (sorted_bam))
    #############################
    print ("Case 1 trim took ", time.time() - start_time, "to run")


    # Take the count of the VSGs that have a mosaic in read 1 (will be the count of the aligned in trimmer that are not 221)
    start_time = time.time()
    
    VSG_count = final_sam + "VSG_count_1.txt"
    os.system("htseq-count -m union -f bam --nonunique all --stranded=no %s %s > %s" % (case1, gff_file, VSG_count))
    
    with open(VSG_count) as f: #Print reads to file
        reader = csv.reader(f, delimiter="\t")
        case1_count = list(reader) #list of the HTSeq count outputs
    print ("Count 1 ", time.time() - start_time, "to run")
    
    
    start_time = time.time()
    case1_mosaic = locations_trim(case1) #arrays of the aligned sequences and positions
    print ("Case 1 locations ", time.time() - start_time, "to run")
    
    #output the aligned file to bam to be used for pysam (this is used for seeing what reads to be considered the later steps)
    sam_file = aligned_read1
    bam_file = sam_file.split(".")[0] + ".bam"
    sorted_bam = bam_file.split(".")[0] + ".sorted.bam"
    remove_list.append(bam_file)
    remove_list.append(sorted_bam)
    remove_list.append(sorted_bam + ".bai")
    
    os.system("samtools view -bS -o %s %s" % (bam_file, sam_file))
    os.system("samtools sort %s -o %s" % (bam_file, sorted_bam))
    os.system("samtools index %s" % (sorted_bam))
    
    samfile = pysam.AlignmentFile(sorted_bam, "rb")
    align_221 = []
    for read in samfile.fetch(target_gene): # find all the reads in the read 1 that align to 221
        align_221.append(str(read.query_name[:-1])) ######## LOOK AT WHOLE READ NAMES
        
        
    
    output_name = final_sam + "case_2.fq" #name of file to be analyzed for case 2
    output = open(output_name, 'w')
    read1 = open(reads_R1, 'rU')
    
    #loop through the reads aligned to 221 and find their corresponding sequence in Read 2
    
    start_time = time.time()
    i = 0
    for title, seq, qual in FastqGeneralIterator(read1):
        read_name = title[:-1] # if the barcode matched a barcode that matched to 221
        
        if read_name in align_221:
            output.write('@' + title + '\n')
            output.write(seq + '\n')
            output.write('+' + '\n')
            output.write(qual + '\n')
            i = i+1
    #print ("New version took", time.time() - start_time, "to run")
    
    output.close()
    
    #find all the reads that align to other VSGs in read 2
    aligned_read2 = final_sam + "aligned_read2.sam"
    unaligned_read2 = final_sam + "unaligned_read2.fq"
    os.system(bowtie_code % (unaligned_read2, reference, output_name, aligned_read2))
    
    case2_locations = locations_gap(aligned_read2, aligned_read1)
    
    start_time = time.time()
    VSG_count = final_sam + "VSG_count_2.txt"
    # Takes the count of the case 2
    os.system("htseq-count -m union -f bam --nonunique all --stranded=no %s %s > %s" % (aligned_read2, gff_file, VSG_count))
    
    with open(VSG_count) as f:
        reader = csv.reader(f, delimiter="\t")
        case2_count = list(reader) #list of the HTSeq count outputs for case 2
    print ("Count 2", time.time() - start_time, "to run")
    
    
    start_time = time.time()
    
    print ("Count 221", time.time() - start_time, "to run")
    # Trims read 2
    start_time = time.time()
    print ("Trimming case 3")
    case3 = trimmer(unaligned_read2)
    
    ###########################  SAVE ALL TRIMMED READ FILES AS BAM
    sam_file = case3
    bam_file = final_sam + "Case3_aligned_trim.bam"
    sorted_bam = bam_file.split(".")[0] + ".sorted.bam"
    
    os.system("samtools view -bS -o %s %s" % (bam_file, sam_file))
    os.system("samtools sort %s -o %s" % (bam_file, sorted_bam))
    os.system("samtools index %s" % (sorted_bam))
    #############################
    
    print ("Case 3 trim ", time.time() - start_time, "to run")
    
    
    print ("Locations case 3")
    start_time = time.time()
    case3_mosaic = locations_trim(case3)
    print ("Case 3 locations", time.time() - start_time, "to run")
    
    VSG_count = final_sam + "VSG_count_3.txt"
    # takes the count of Case 3 (mosaic in read 2)
    os.system("htseq-count -m union -f bam --nonunique all --stranded=no %s %s > %s" % (case3, gff_file, VSG_count))
    with open(VSG_count) as f:
        reader = csv.reader(f, delimiter="\t")
        case3_count = list(reader) #list of the HTSeq count outputs

    # For each of the other VSGs, the frequency of mosaics is the addition of the count for the three cases
    mosaic_count = {}
    
    for VSG_gene in reference_VSG:
        i = 0
        rows = len(case3_count)
        while i < rows:
            if VSG_gene == case3_count[i][0]:
                if VSG_gene == target_gene:
                    mosaic_count[VSG_gene] = int(case2_count[i][1])
                else:
                    mosaic_count[VSG_gene] = int(case1_count[i][1]) + int(case2_count[i][1]) + int(case3_count[i][1])
            i = i + 1
    
    #print(mosaic_count)

    remove_list.append(final_sam + "_output_readnames.fq")
    remove_list.append(final_sam + "case_2.fq")
    remove_list.append(final_sam + "trimmed_unaligned.fq")
    remove_list.append(final_sam + "unaligned_read1.fq")
    remove_list.append(final_sam + "unaligned_read2.fq")
    remove_list.append(final_sam + "unaligned.fq")
    
    remove_list.append(final_sam + "aligned_read1.sam")
    remove_list.append(final_sam + "aligned_read2.sam")
    remove_list.append(final_sam + "aligned_100.sam")
    
    ##################################################

    final_output = open(final_sam + "output.txt" , 'w')
    location_output = open(final_sam + "output_locations.txt" , 'w')
    name_output = open(final_sam + "output_readname.txt" , 'w')
    
    final_output.write(bowtie_code + '\n')
    
    final_output.write("Total reads: " + str(validate_count) + '\n')
    
    sum_count = 0
    for entry in mosaic_count:
        sum_count = sum_count + int(mosaic_count[entry])
    
    final_output.write("Total aligned: " + str(sum_count) + '\n')
    #final_output.write("Unaligned estimate: " + str(validate_count["Total trim 2 unaligned"] + validate_count["Total trim 1 unaligned"]) + '\n')
    
    #print(mosaic_count)
    #print(mosaic_count[target_gene])
    
    if sum_count > 0:
        final_output.write("Nonmosaic %s" % (target_gene) + "\t" + str(mosaic_count[target_gene]) + '\n') # print the frequency of non mosaics
    
    # print to the file and then to the screen
    print("\nFinal output written to %s output.txt" % (final_sam))
    print("\n%s") % (target_gene)
    #print(mosaic_count[target_gene])
    i = 0
    j = 0
    
    for entry in mosaic_count:
        if entry != target_gene:
            final_output.write(entry + '\t' + str(mosaic_count[entry]) + '\n')
    
    # write the mosaics and positions
    #final_output.write("\ncase 1")
    
    pos_dict = {}
    
    
    i = 0
    
    for reads in case1_mosaic:
        try:
            a = case1_mosaic[reads][1]
        except:
            pass
        else:
            name_output.write('@' + reads)
            name_output.write('\t' + case1_mosaic[reads][1])
            name_output.write('\t' + case1_mosaic[reads][0])
            name_output.write('\t' + case1_mosaic[reads][2] + '\n')
            try:
                a = pos_dict[case1_mosaic[reads][0]]
            except:
                pos_dict[case1_mosaic[reads][0]] = {}
                pos_dict[case1_mosaic[reads][0]][case1_mosaic[reads][1]] = 1
            else:
                try:
                    a = pos_dict[case1_mosaic[reads][0]][case1_mosaic[reads][1]]
                except:
                    pos_dict[case1_mosaic[reads][0]][case1_mosaic[reads][1]] = 1
                else:
                    pos_dict[case1_mosaic[reads][0]][case1_mosaic[reads][1]] = pos_dict[case1_mosaic[reads][0]][case1_mosaic[reads][1]] + 1
        i = i + 1
    
    #name_output.write("\n\ncase 2")
    i = 0
    for reads in case2_locations:
        #final_output.write('\n' + case2_locations[reads][0])
        #final_output.write('\t' + case2_locations[reads][1])
        try:
            a = pos_dict[case2_locations[reads][1]][case2_locations[reads][0]]
        except:
            pos_dict[case2_locations[reads][1]] = {}
            pos_dict[case2_locations[reads][1]][case2_locations[reads][0]] = 1
        else:
            name_output.write('@' + reads)
            name_output.write('\t' + case2_locations[reads][0])
            name_output.write('\t' + case2_locations[reads][1] + '\n')
            try:
                a = pos_dict[case2_locations[reads][1]][case2_locations[reads][0]]
            except:
                pos_dict[case2_locations[reads][1]][case2_locations[reads][0]] = 1
            else:    
                pos_dict[case2_locations[reads][1]][case2_locations[reads][0]] = pos_dict[case2_locations[reads][1]][case2_locations[reads][0]] + 1
        i = i + 1
    
    #final_output.write("\n\ncase 3")
    i = 0
    for reads in case3_mosaic:
        try:
            a = case3_mosaic[reads][1]
        except:
            pass
        else:
            name_output.write('@' + reads)
            name_output.write('\t' + case3_mosaic[reads][1])
            name_output.write('\t' + case3_mosaic[reads][0])
            name_output.write('\t' + case3_mosaic[reads][2] + '\n')
            try:
                a = pos_dict[case3_mosaic[reads][0]][case3_mosaic[reads][1]]
            except:
                pos_dict[case3_mosaic[reads][0]] = {}
                pos_dict[case3_mosaic[reads][0]][case3_mosaic[reads][1]] = 1
            else:
                try:
                    a = pos_dict[case3_mosaic[reads][0]][case3_mosaic[reads][1]]
                except:
                    pos_dict[case3_mosaic[reads][0]][case3_mosaic[reads][1]] = 1
                else:
                    pos_dict[case3_mosaic[reads][0]][case3_mosaic[reads][1]] = pos_dict[case3_mosaic[reads][0]][case3_mosaic[reads][1]] + 1
        i = i + 1
    
    
    location_output.write("Location")
    for VSG in reference_VSG:
        if VSG != target_gene:
            location_output.write('\t' + VSG)
        
    
    
    pos_list = sorted(pos_dict)
    #print pos_list
    
    for position in pos_list:
        #print position
        location_output.write('\n' + position)
        for VSG in reference_VSG:
            if VSG != target_gene:
                try:
                    a = pos_dict[position][VSG]
                except:
                    location_output.write('\t' + '0')
                else:
                    location_output.write('\t' + str(a))
                
    #print(remove_list)
    remove_files(remove_list)

def demultiplex(bc, R1, R2):
    barcode_file = open(bc, 'rU')
    barcode_dict = {}
    
    for line in barcode_file:
        sample_data = line.strip().split('\t')
        barcode_dict[sample_data[1]+"+"+sample_data[2]] = sample_data[0]
    
    len_dict = []
    for barcode in barcode_dict:
        for read, fastq in [("R1", R1), ("R2",R2)]:
            outfile = open(barcode_dict[barcode]+"_"+str(barcode)+"_"+read+'.fq', 'w')
            with open(fastq) as in_handle:
                for title, seq, qual in FastqGeneralIterator(in_handle):
                    if title.split(":")[-1][:-8] == barcode:
                        UMI = title[-8:]
                        #if len(title[:-8]) == 68 :
                        len_dict.append(len(title[:-8]))
                        
                        readname = title.split(" ")[0]+'_'+UMI+"\\"+read[1] 
                        #print readname
                        #print len(readname)
                        #print "@%s\n%s\n+\n%s\n" % (readname, seq, qual)
                        outfile.write("@%s\n%s\n+\n%s\n" % (readname, seq, qual))

def names(barcodes):
    barcode_file = open(barcodes, "r")
    
    demultiplex_out = {}
    consolidate_out = {}
    
    for line in barcode_file:
        #print line
        base_name = line.split("\t")[0] + "_" + line.split("\t")[1] + "+" + line.split("\t")[2]
        base_name = base_name.strip()
        #print base_name + "tab"
        
        demultiplex_out[base_name] = []
        demultiplex_out[base_name].append(base_name + "_R1.fq")
        demultiplex_out[base_name].append(base_name + "_R2.fq")
        
        consolidate_out[base_name] = []
        consolidate_out[base_name].append(base_name + "_consolidated_R1.fq")
        consolidate_out[base_name].append(base_name + "_consolidated_R2.fq")
    
    return demultiplex_out, consolidate_out

def consolidate(read1, read2, Run_name, min_qual, min_freq):
    start_program = time.time()
    logger = logging.getLogger('root')
    
    def assure_match(molecular_id_dict, junk_file):
        junk = open(junk_file, 'w') #Creat junk file to catch all the reads that are not consolidated
        junk_sequence = junk_file[:-5] + "_junk_sequence"
        junk_reads = open(junk_sequence, 'w')
        
        junk.write("%s\t%s\t%s\n" % ("Molecular id", "Read Start", "Reads"))
        junk_dict = []
        for molecular_id in molecular_id_dict:
            for read_start in molecular_id_dict[molecular_id].keys(): #Iterate through the Readstarts within each UMI
                if len(molecular_id_dict[molecular_id][read_start]) > 1: #If the group is larger than just one read with a UMI and readstart
                    ratio = []
                    matching = {}
                    length = len(molecular_id_dict[molecular_id][read_start]) #size of the group
                    number = 0
                    
                    # matching[number] = {}
                    # matching[number] = [molecular_id_dict[molecular_id][read_start][0]]
                    
                    while number < length: # go through the whole group starting at the first read
                        matching[number] = {}
                        matching[number] = [molecular_id_dict[molecular_id][read_start][number]] #add the first read to the matching dictionary under new integer index
                        #print number
                        i = 0
                        delete = [] #List of reads to be deleted as they are similar to the other reads
                        while i < length: #compare the current read to all others and if any is similar above the 98% threshhold then add to the delete file
                            if i != number: #so as to not compare a read to itself
                                read_length = min(len(molecular_id_dict[molecular_id][read_start][number].seq), len(molecular_id_dict[molecular_id][read_start][i].seq)) #compare only the shared bases
                                #print molecular_id_dict[molecular_id][read_start][number].seq[0:read_length]
                                #print molecular_id_dict[molecular_id][read_start][i].seq[0:read_length]
                                
                                s = SequenceMatcher(None, molecular_id_dict[molecular_id][read_start][number].seq[0:read_length], molecular_id_dict[molecular_id][read_start][i].seq[0:read_length]) #compare the two reads
                                #print s.ratio()
                                
                                if s.ratio() >= 0.98: #if the similarity is larger than 98%
                                    delete.append(i) #add the compared read to delete
                                    if not number in matching: 
                                        matching[number] = {}
                                        matching[number] = [molecular_id_dict[molecular_id][read_start][i]]
                                    else:
                                        matching[number].append(molecular_id_dict[molecular_id][read_start][i]) #add the read to the group
                                elif s.ratio() < 0.98 and ((molecular_id + read_start) not in junk_dict): #so as to only add the not matching reads once
                                    junk.write("%s\t%s\t" % (molecular_id, read_start))
                                    for read in molecular_id_dict[molecular_id][read_start]:
                                        junk.write(read.name+';')
                                    junk.write('\n')
                                    junk_dict.append(molecular_id + read_start) # make a list of the reads that are not in the group
                                    
                                    junk_reads.write(molecular_id + "\t" + read_start + "\n")
                                    
                                    for read in molecular_id_dict[molecular_id][read_start]:
                                        junk_reads.write(read.name+'; '+ read.seq + "\n")
                                
                                ratio.append(s.ratio())
                                print("")
                            i = i + 1
                        number = number + 1
                        #print delete
                        
                        for read in sorted(delete, reverse=True):
                            del(molecular_id_dict[molecular_id][read_start][read]) #remove the matching reads from the list
                        
                        length = len(molecular_id_dict[molecular_id][read_start])   #changes the length of the list so you can move to the next read
                    
                    del(molecular_id_dict[molecular_id][read_start]) #remove the group of reads
                    
                    for number in matching: #each number represents a different group of matched reads
                        i = 0
                        for read in matching[number]:
                            new_readname = read_start + str(number) #add an index to the end of the readstart
                            
                            if not molecular_id in molecular_id_dict: #add the reads back to the dictionary with the slightly new groups (these will fall under the same UMI but will have an edited read start)
                                molecular_id_dict[molecular_id] = {}
                                molecular_id_dict[molecular_id][new_readname]= [matching[number][i]]
                            else:
                                if not new_readname in molecular_id_dict[molecular_id]:
                                    molecular_id_dict[molecular_id][new_readname] = [matching[number][i]]
                                else:
                                    molecular_id_dict[molecular_id][new_readname].append(matching[number][i])
                            i = i + 1
    
    def sort_reads_twosided(fastq_file_1, fastq_file_2, junk_file):
        infile = HTSeq.FastqReader(fastq_file_1) #Looking at read 1 first (far end, this is where we are taking the first 8 base pairs)
        read_num = 0
    
        molecular_id_dict = {}
    
        for read in infile: #loop through the reads
            read_name = read.name.split("_")[0] 
            #print read_name
            sample_id = str(fastq_file_1).split("_")[-2] 
            #print sample_id
            molecular_id = read.name.split("_")[1][0:8] #unique UMI
            read_start = read.seq[0:8] #The first 8 base pairs of Read 1
    
            if not molecular_id in molecular_id_dict: #adding the molecular id to a dictionary (UMI:[Read start1:[Reads], Readstart2:[Reads], Readstart3:[Reads]])
                molecular_id_dict[molecular_id] = {} #create new index for UMI
                molecular_id_dict[molecular_id][read_start]= [read]
            else:
                if not read_start in molecular_id_dict[molecular_id]: #create mew index for read start
                    molecular_id_dict[molecular_id][read_start] = [read]
                else:
                    molecular_id_dict[molecular_id][read_start].append(read) #Add read to existing UMI and readstart
        #print molecular_id_dict["GTAGTTAT"]
        #print molecular_id_dict["GTAGTTAT"]["ATCCGAGC"][0].name
        
        assure_match(molecular_id_dict, junk_file)
        
        read2 = HTSeq.FastqReader(fastq_file_2)
        molecular_id_dict_2 = {}
        
        for read in read2:
            read_name = read.name.split("\\")[0] # this should match the matching read in file1
            #print read_name
            molecular_id = read.name.split("_")[1][0:8]
            read_start_actual = ""
            for read_start in molecular_id_dict[molecular_id]:
                i = 0
                length = len(molecular_id_dict[molecular_id][read_start])
                while i < length:
                    if read_name == molecular_id_dict[molecular_id][read_start][i].name.split("\\")[0]: #loop through the groups in read1 and find the matching reads in read2
                        read_start_actual = read_start
                    i = i + 1
            
            if read_start_actual == "":
                print("Error: not even reads")
                
            if not molecular_id in molecular_id_dict_2: #add these reads to a new dictionary
                molecular_id_dict_2[molecular_id] = {}
                molecular_id_dict_2[molecular_id][read_start_actual]= [read]
            else:
                if not read_start_actual in molecular_id_dict_2[molecular_id]:
                    molecular_id_dict_2[molecular_id][read_start_actual] = [read]
                else:
                    molecular_id_dict_2[molecular_id][read_start_actual].append(read)
        
        #print molecular_id_dict["GTAGTTAT"]
        #print molecular_id_dict_2["GTAGTTAT"]
        return molecular_id_dict, molecular_id_dict_2
    
    def sort_reads(fastq_file):
        infile = HTSeq.FastqReader(fastq_file)
        read_num = 0
        #bin_name = [] # Will consist of the molecular_id followed by the incoming 2nd header field (which includes sample_id)
        #bin_reads = []
        molecular_id_dict = {}
        #import itertools
        #for read in itertools.islice( infile, 10000 ): # Read the first few reads for testing
        for read in infile:
            read_name = read.name.split("_")[0]
            #print read_name
            sample_id = str(fastq_file).split("_")[-2]
            #print sample_id
            molecular_id = read.name.split("_")[1][0:8]
            read_start = read.seq[0:8]
            #print molecular_id +"molecular_id"
            #print read_start +"read_start"
            if not molecular_id in molecular_id_dict:
                molecular_id_dict[molecular_id] = {}
                molecular_id_dict[molecular_id][read_start]= [read]
            else:
                if not read_start in molecular_id_dict[molecular_id]:
                    molecular_id_dict[molecular_id][read_start] = [read]
                else:
                    molecular_id_dict[molecular_id][read_start].append(read)
        #print molecular_id_dict["GTAGTTAT"]
        #print molecular_id_dict["GTAGTTAT"]["GGCTTGTG"]
        return molecular_id_dict
    
    def consolidate_position(bases, quals, min_qual, min_freq):
        num = {}
        qual = {}
        num['A'] = num['C'] = num['G'] = num['T'] = num['N'] = 0
        qual['A'] = qual['C'] = qual['G'] = qual['T'] = qual['N'] = 0
        for bb, qq in zip(bases, quals):
            if qq > min_qual:
                num[bb] += 1
            if qq > qual[bb]:
                qual[bb] = qq
        most_common_base = max(num.iterkeys(), key=(lambda key: num[key]))
        freq = float(num[most_common_base]) / len(bases)
        if freq > min_freq:
            return True, most_common_base, qual[most_common_base]
        else:
            return False,'N', 0
        
    def consolidate(fastq_file, consolidated_fastq_file, read_data_file, min_qual, min_freq, bins, suffix):
        outfile = open(consolidated_fastq_file, 'w')
        readdatafile = open(read_data_file + suffix, 'w')
        readdatafile.write('Molecular_id\tnumber\tnumber_reads_used\treads\tRead_start\n')
        #bins = read_bins(fastq_file)
        #bins = sort_reads(fastq_file)
        #next(bins)
    
        #outfile.write("SHANE")
    
        num_input_reads = 0
        num_consolidated_reads = 0
        num_successes = 0 # Bases with successful consolidation
        num_bases = 0
        for id in bins.iterkeys():
            num_start = 0
            for start in bins[id].iterkeys():
                reads = bins[id][start]
                num_start += 1
                num_input_reads += len(reads)
                num_consolidated_reads += 1
                # Get all the bases and quals in the read
                read_bases = zip(*[list(read.seq) for read in reads])
                read_quals = zip(*[list(read.qual) for read in reads])
                # Iterate position by position
                consolidation_sucess, cons_seq, cons_qual = zip(*[consolidate_position(bases, quals, min_qual, min_freq) for bases, quals in zip(read_bases, read_quals)])
                # Count consolidation successes and failures
                num_successes += sum(consolidation_sucess)
                num_bases += len(consolidation_sucess)
                # Write consolidated FASTQ read
                outfile.write('@%s_%d_%dreads\%s\n' % (id, num_start,len(reads), suffix)) # Header: Molecular id, position for this molecular id, number of reads for consolidation
                outfile.write(''.join(cons_seq) +'\n')
                outfile.write('+\n')
                outfile.write(''.join([chr(q+33) for q in cons_qual]) + '\n')
                readdatafile.write('%s\t%d\t%d\t%s\t' % (id, num_start,len(reads), start))
                for read in reads:
                    readdatafile.write(read.name+';')
                readdatafile.write('\n')
    
        logger.info("Read %d input reads", num_input_reads)
        logger.info("Wrote %d consolidated reads", num_consolidated_reads)
        logger.info("Successfully consolidated %d bases out of %d (%.2f%%)", num_successes, num_bases, 100*float(num_successes)/num_bases)
        outfile.close()


    fastq_file_1 = read1
    fastq_file_2 = read2
    
    
    consolidated_fastq_file = Run_name + "_" + "consolidated_R1.fq" #Final file names
    consolidated_fastq_file_2 = Run_name + "_" + "consolidated_R2.fq" 
    read_data_file = Run_name + "_read_data"
    
    junk = Run_name + "_junk" #junk file

    min_qual = int(min_qual) #the minimum quality for consolidation (15)
    min_freq = float(min_freq) #minimum frequency of base pair (0.9)

    molecular_id_dict, molecular_id_dict_2 = sort_reads_twosided(fastq_file_1, fastq_file_2, junk)
    
    print ("Consolidate time", time.time() - start_program)
    
    consolidate(fastq_file_1, consolidated_fastq_file, read_data_file, min_qual, min_freq, molecular_id_dict, "1") #consolidate each of the two files
    consolidate(fastq_file_2, consolidated_fastq_file_2, read_data_file, min_qual, min_freq, molecular_id_dict_2, "2")
    
    print ("Final time", time.time() - start_program)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-bc',help='tab-delimited text file with sample name, i7 barcode sequence, and i5 barcode sequence ', action="store", dest="bc", default='')
    parser.add_argument('-R1',help='R1 fastq file', action="store", dest="R1", default='')
    parser.add_argument('-R2',help='R2 fastq file', action="store", dest="R2", default='')
    parser.add_argument('-min_qual',help='Minimum quality for consolidation. Default = 15', action="store", dest="min_qual", default='15')
    parser.add_argument('-min_freq',help='Minimum frequency for consolidation. Default = 0.9', action="store", dest="min_freq", default='0.9')
    parser.add_argument('-ref',help='Reference name for the genome (.ebwt file name)', action="store", dest="reference", default='')
    parser.add_argument('-genome',help='The genome in fasta file format (.fasta file)', action="store", dest="genome_name", default='')
    parser.add_argument('-gff',help='VSG informational file. (.gff file)', action="store", dest="gff_file", default='')
    parser.add_argument('-target',help='The name of the target gene. (ex. VSG_221 or Tb427VSG-2)', action="store", dest="target_gene", default='')
    parser.add_argument('--extra',help='Default is to remove the extra files. These files are the VSG_count for the three cases, Val_1 and Val_2 for the demultiplexing and trim.bam files --extra indicates keeping them', dest="extra", action='store_true')
    parser.add_argument('--no-extra',help='Default is to remove the extra files. --no-extra indicates removing extra files', dest="extra", action='store_false')
    parser.set_defaults(extra=False)
    parser.add_argument('--junk',help='Default is to remove the junk files. --junk indicates keeping them', dest="junk", action='store_true')
    parser.add_argument('--no-junk',help='Default is to remove the junk files. --no-junk indicates removing junk files', dest="junk", action='store_false')
    parser.set_defaults(junk=False)
    
    
    arguments = parser.parse_args()
    
    # inputfile_R1 = str(argv[1]) # unzipped file
    # inputfile_R2 = str(argv[2]) # unzipped file read 2
    # barcodes = str(argv[3])
    # min_qual = str(argv[4])
    # min_freq = str(argv[5])
    #reference = str(argv[6])
    #genome_name = str(argv[7])
    #gff_file = str(argv[8])
    #target_gene = str(argv[9])
    
    print("TRIMMING")
    
    os.system("trim_galore --dont_gzip --paired --adapter ACTCCAGTCAC --adapter2 AGATCGGAAGAGC --trim1 %s %s" % (arguments.R1, arguments.R2))
    
    trim_out_R1 = arguments.R1.split(".")[0] + "_val_1.fq"
    trim_out_R2 = arguments.R2.split(".")[0] + "_val_2.fq"
    
    print("DEMULTIPLEXING")
    
    demultiplex(arguments.bc, trim_out_R1, trim_out_R2)
    
    demultiplex_names, consolidate_names = names(arguments.bc)
    
    print(consolidate_names)
    
    print("CONSOLIDATING")
    
    for name in demultiplex_names:
        #print demultiplex_names[name][0]
        consolidate(demultiplex_names[name][0], demultiplex_names[name][1], name, arguments.min_qual, arguments.min_freq)

    print("VSG_seq")

    for name in consolidate_names:
        #print(name)
        VSG(consolidate_names[name][0], consolidate_names[name][1], arguments.reference, name, arguments.genome_name, arguments.gff_file, arguments.target_gene)
    
    if not arguments.extra:
        print("Removing extra Files")
        os.remove(trim_out_R1)
        os.remove(trim_out_R2)
        for name in consolidate_names:
            base = name.split("_consolidate")[0]
            count1 = base + "VSG_count_1.txt"
            count2 = base + "VSG_count_2.txt"
            count3 = base + "VSG_count_3.txt"
            case1 = base + "Case1_aligned_trim.bam"
            case3 = base + "Case3_aligned_trim.bam"
            os.remove(count1)
            os.remove(count2)
            os.remove(count3)
            os.remove(case1)
            os.remove(case3)
    
    if not arguments.junk:
        print("Removing Junk Files")
        for name in consolidate_names:
            base = name.split("_consolidate")[0]
            trim_junk = base + "_junk.fq"
            consolidate_junk = base + "_junk"
            consolidate_seq = base + "_junk_sequence"
            os.remove(trim_junk)
            os.remove(consolidate_junk)
            os.remove(consolidate_seq)

if __name__ == '__main__':
    main()
