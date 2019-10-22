###### VSG AMP-Seq Counting
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
logger = logging.getLogger('root')

from Bio import pairwise2 as pw2

from difflib import SequenceMatcher



def locations_trim(reads, case): #input the sam file from the trimmer
    sam = reads
    bam = sam.split(".")[0] + ".bam"
    sorted_bam = bam.split(".")[0] + ".sorted.bam"
    
    remove_list.append(bam)
    remove_list.append(sorted_bam)
    remove_list.append(sorted_bam + ".bai")
    
    os.system("samtools view -bS -o %s %s" % (bam, sam)) # Turn the sam file into a sorted BAM with an index
    os.system("samtools sort %s -o %s" % (bam, sorted_bam))
    os.system("samtools index %s" % (sorted_bam))
    
    try:
        bamfile = pysam.AlignmentFile(sorted_bam, "rb") #uses a bam file with an index
    except:
        print("error in reading bam files")
        global error
        global error_case1
        global error_case2
        global error_case3
        error = True
        if case == "case1":
            error_case1 = True
        elif case == "case2":
            error_case2 = True
        elif case == "case3":
            error_case3 = True
        return
    
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

def locations_gap(reads, case2, case): # input the sam file and the case 2 fq file
    sam = reads
    bam = sam.split(".")[0] + ".bam"
    sorted_bam = bam.split(".")[0] + ".sorted.bam"
    
    remove_list.append(bam)
    remove_list.append(sorted_bam)
    remove_list.append(sorted_bam + ".bai")
    
    os.system("samtools view -bS -o %s %s" % (bam, sam)) # Turn the sam file into a sorted BAM with an index
    os.system("samtools sort %s -o %s" % (bam, sorted_bam))
    os.system("samtools index %s" % (sorted_bam))
    
    try:
        bamfile = pysam.AlignmentFile(sorted_bam, "rb") #uses a bam file with an index
    except:
        print("error in reading bam file")
        global error
        global error_case1
        global error_case2
        global error_case3
        error = True
        if case == "case1":
            error_case1 = True
        elif case == "case2":
            error_case2 = True
        elif case == "case3":
            error_case3 = True
        return
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

def count_reads(filename):
    i = 0
    for title in FastqGeneralIterator(filename):
        i = i+1
    validate_count = i
    return validate_count

##########################################################################
    print("starting " + final_sam + "........")

def read_count(case_name, gff_file, output_name, case):
    failure = os.system("python -m HTSeq.scripts.count -m union -f bam --nonunique all --stranded=no %s %s > %s" % (case_name, gff_file, output_name))
    print(failure)
    if failure != 0:
        print("error")
        global error
        global error_case1
        global error_case2
        global error_case3
        error = True
        if case == "case1":
            error_case1 = True
        elif case == "case2":
            error_case2 = True
        elif case == "case3":
            error_case3 = True
    else:
        with open(output_name) as f: #Print reads to file
            reader = csv.reader(f, delimiter="\t")
            case_count = list(reader) #list of the HTSeq count outputs
    
        return case_count

def names(barcodes):
    # Takes an input of the barcode file used in the experiment and makes a list of files for each of the different experimental
    # conditions for the demultiplexed and then consolidated files
    
    barcode_file = open(barcodes, "r")
    
    file_names = []
    
    for line in barcode_file:
        #print line
        base_name = line.split("\t")[0] + "_" + line.split("\t")[1] + "+" + line.split("\t")[2]
        base_name = base_name.strip()
        #print base_name + "tab"
        
        file_names.append(base_name)
    
    return file_names

def error_analyze(error_case1, error_case2, error_case3, final_sam):
    #inputs
    case1 = final_sam + "aligned_100case1.sam" # From trimmer
    case3 = final_sam + "aligned_100case3.sam" # From trimmer
    aligned_read1 = final_sam + "aligned_read1.sam" # Bowtie alignment of read1
    aligned_read2 = final_sam + "aligned_read2.sam" # Bowtie alignment of read2
    #outputs
    VSG_count1 = final_sam + "VSG_count_1.txt"
    VSG_count2 = final_sam + "VSG_count_2.txt"
    VSG_count3 = final_sam + "VSG_count_3.txt"
    
    if error_case1 == False and error_case2 == True and error_case3 == True: #Case of only Case 1 being without errors (this would be only the read with the primer being funcitonal)
        case1_count = read_count(case1, gff_file, VSG_count1, "case1")
        case1_mosaic = locations_trim(case1, "case1") #arrays of the aligned sequences and positions
        mosaic_count = {}
        for VSG_gene in reference_VSG:
            i = 0
            rows = len(case1_count)
            while i < rows:
                if VSG_gene == case1_count[i][0]:
                    if VSG_gene == target_gene:
                        print("Full VSG")
                        #mosaic_count[VSG_gene] = int(case2_count[i][1])
                    else:
                        mosaic_count[VSG_gene] = int(case1_count[i][1])
                i = i + 1
        
        final_output = open(final_sam + "output.txt" , 'w')
        location_output = open(final_sam + "output_locations.txt" , 'w')
        name_output = open(final_sam + "output_readname.txt" , 'w')
        
        final_output.write(bowtie_code + '\n')
        
        sum_count = 0
        for entry in mosaic_count:
            sum_count = sum_count + int(mosaic_count[entry])
        
        final_output.write("Total aligned: " + str(sum_count) + '\n')
        
        try:
            final_output.write("Nonmosaic %s" % (target_gene) + "\t" + str(mosaic_count[target_gene]) + '\n') # print the frequency of non mosaics
        except:
            final_output.write("Nonmosaic %s" % (target_gene) + "\t" + '0' + '\n') # If no mosaics are present then just print 0
        
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
    elif error_case1 == True and error_case2 == True and error_case3 == True: # All three cases with errors
        print("Too many errors")
    else:
        if error_case1 == False:
            case1_count = read_count(case1, gff_file, VSG_count1, "case1")
            case1_mosaic = locations_trim(case1, "case1") #arrays of the aligned sequences and positions
        if error_case2 == False:
            case2_locations = locations_gap(aligned_read2, aligned_read1, "case2")
            case2_count = read_count(aligned_read2, gff_file, VSG_count1, "case2")
        if error_case3 == False:
            case3_mosaic = locations_trim(case3, "case3")
            case3_count = read_count(case3, gff_file, VSG_count1, "case3")

        # For each of the other VSGs, the frequency of mosaics is the addition of the count for the three cases
        mosaic_count = {}
        
        if error_case1 == False:
           for VSG_gene in reference_VSG:
                i = 0
                rows = len(case1_count)
                while i < rows:
                    if VSG_gene == case1_count[i][0]:
                        if VSG_gene == target_gene:
                            print("nope")
                            #mosaic_count[VSG_gene] = int(case2_count[i][1])
                        else:
                            mosaic_count[VSG_gene] = int(case1_count[i][1])
                    i = i + 1
        if error_case2 == False:
            for VSG_gene in reference_VSG:
                i = 0
                rows = len(case2_count)
                while i < rows:
                    if VSG_gene == case2_count[i][0]:
                        if VSG_gene == target_gene:
                            mosaic_count[VSG_gene] = int(case2_count[i][1])
                        else:
                            try:
                                mosaic_count[VSG_gene] = mosaic_count[VSG_gene] + int(case2_count[i][1])
                            except:
                                mosaic_count[VSG_gene] = int(case2_count[i][1])
                                
                    i = i + 1
        if error_case3 == False:
            for VSG_gene in reference_VSG:
                i = 0
                rows = len(case3_count)
                while i < rows:
                    if VSG_gene == case2_count[i][0]:
                        if VSG_gene == target_gene:
                            print("nope")
                        else:
                            try:
                                mosaic_count[VSG_gene] = mosaic_count[VSG_gene] + int(case3_count[i][1])
                            except:
                                mosaic_count[VSG_gene] = int(case3_count[i][1])
        
        final_output = open(final_sam + "output.txt" , 'w')
        location_output = open(final_sam + "output_locations.txt" , 'w')
        name_output = open(final_sam + "output_readname.txt" , 'w')
        
        final_output.write(bowtie_code + '\n')
        
        #final_output.write("Total reads: " + str(validate_count) + '\n')
        
        sum_count = 0
        for entry in mosaic_count:
            sum_count = sum_count + int(mosaic_count[entry])
        
        final_output.write("Total aligned: " + str(sum_count) + '\n')
    
        try:
            final_output.write("Nonmosaic %s" % (target_gene) + "\t" + str(mosaic_count[target_gene]) + '\n') # print the frequency of non mosaics
        except:
            final_output.write("Nonmosaic %s" % (target_gene) + "\t" + '0' + '\n') # If no mosaics are present then just print 0
        
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
        if error_case1 == False:
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
        if error_case2 == False:
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
        if error_case3 == False:
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
        #remove_files(remove_list)
        print("werid error")

def analyze_reads(final_sam):
    #inputs
    case1 = final_sam + "aligned_100case1.sam" # From trimmer
    case3 = final_sam + "aligned_100case3.sam" # From trimmer
    aligned_read1 = final_sam + "aligned_read1.sam" # Bowtie alignment of read1
    aligned_read2 = final_sam + "aligned_read2.sam" # Bowtie alignment of read2
    #outputs
    VSG_count1 = final_sam + "VSG_count_1.txt"
    VSG_count2 = final_sam + "VSG_count_2.txt"
    VSG_count3 = final_sam + "VSG_count_3.txt"
    
    global error_case1
    global error_case2
    global error_case3
    error_case1 = False
    error_case2 = False
    error_case3 = False
    
    print("case1")
    case1_count = read_count(case1, gff_file, VSG_count1, "case1")
    case1_mosaic = locations_trim(case1, "case1") #arrays of the aligned sequences and positions
    print("case2")
    case2_locations = locations_gap(aligned_read2, aligned_read1, "case2")
    case2_count = read_count(aligned_read2, gff_file, VSG_count2, "case2")
    print("case3")
    case3_mosaic = locations_trim(case3, "case3")
    case3_count = read_count(case3, gff_file, VSG_count3, "case3")
    
    if error == True:
        print("There has been an error in analyzing %s" % (final_sam))
        if error_case1 == True:
            print("error_case1")
        if error_case2 == True:
            print("error_case2")
        if error_case3 == True:
            print("error_case3")
        error_analyze(error_case1, error_case2, error_case3, final_sam)
        return
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
    
    final_output = open(final_sam + "output.txt" , 'w')
    location_output = open(final_sam + "output_locations.txt" , 'w')
    name_output = open(final_sam + "output_readname.txt" , 'w')
    
    final_output.write(bowtie_code + '\n')
    
    #final_output.write("Total reads: " + str(validate_count) + '\n')
    
    sum_count = 0
    for entry in mosaic_count:
        sum_count = sum_count + int(mosaic_count[entry])
    
    final_output.write("Total aligned: " + str(sum_count) + '\n')

    try:
        final_output.write("Nonmosaic %s" % (target_gene) + "\t" + str(mosaic_count[target_gene]) + '\n') # print the frequency of non mosaics
    except:
        final_output.write("Nonmosaic %s" % (target_gene) + "\t" + '0' + '\n') # If no mosaics are present then just print 0
    
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
    #remove_files(remove_list)
def remove_files(remove):
    print ("Removing files .......")
    remove = list(set(remove))
    for x in remove:
        #print x
        try:
            os.remove(x)
        except:
            print("Already gone")

def main():
    start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('-name',help='tab-delimited text file with sample name, i7 barcode sequence, and i5 barcode sequence ', action="store", dest="name", default='')
    parser.add_argument('-gff',help='VSG informational file. (.gff file)', action="store", dest="gff_file", default='')
    parser.add_argument('-target',help='The name of the target gene. (ex. VSG_221 or Tb427VSG-2)', action="store", dest="target_gene", default='')
    parser.add_argument('-ref',help='Reference name for the genome (.ebwt file name)', action="store", dest="reference", default='')
    parser.add_argument('-genome',help='The genome in fasta file format (.fasta file)', action="store", dest="genome_name", default='')
    parser.add_argument('-bc',help='tab-delimited text file with sample name, i7 barcode sequence, and i5 barcode sequence ', action="store", dest="bc", default='')
    parser.add_argument('--junk',help='Default is to remove the junk files. --junk indicates keeping them', dest="junk", action='store_true', default = False)
    arguments = parser.parse_args()
    
    global remove_list
    remove_list = []
    global target_gene
    target_gene = arguments.target_gene
    global reference_VSG
    reference_VSG = [] #list of all the VSG genes
    for record in SeqIO.parse(arguments.genome_name, "fasta"):
        reference_VSG.append(record.id)
    global bowtie_code
    bowtie_code = "bowtie -p 20 -v 1 -m 1 -X 1500 --quiet --un %s -S %s %s > %s"
    global gff_file
    gff_file = arguments.gff_file
    global error
    error = False
    global error_case1
    error_case1 = False
    global error_case2
    error_case2 = False
    global error_case3
    error_case3 = False
    
    file_names = names(arguments.bc)
    error_dict = []
    analyzed = []
    
    for name in file_names:
        error = False
        error_case1 = False
        error_case2 = False
        error_case3 = False
        print(name)
        analyze_reads(name)
        if error == True:
            error_dict.append(name)
        else:
            analyzed.append(name)
    
    print("The following reads have been analyzed")
    for read in analyzed:
        print (read)
    print("The following read files contain an error")
    for read in error_dict:
        print (read)
    
    if not arguments.junk:
        print("Removing Junk Files")
        for name in file_names:
            base = name.split("_consolidate")[0]
            trim_junk = base + "_junk.fq"
            consolidate_junk = base + "_junk"
            consolidate_seq = base + "_junk_sequence"
            remove_list.append(trim_junk)
            remove_list.append(consolidate_junk)
            remove_list.append(consolidate_seq)
    
    for name in file_names:
        remove_list.append(name + "Case1_aligned_trim.bam")
        remove_list.append(name + "Case1_aligned_trim.sorted.bam")
        remove_list.append(name + "Case1_aligned_trim.sorted.bam.bai")
        remove_list.append(name + "Case3_aligned_trim.bam")
        remove_list.append(name + "Case3_aligned_trim.sorted.bam")
        remove_list.append(name + "Case3_aligned_trim.sorted.bam.bai")
        
        remove_list.append(name + "trimmed_unaligned.fq")
        remove_list.append(name + "unaligned.fq")
        remove_list.append(name + "unaligned_read1.fq")
        remove_list.append(name + "unaligned_read2.fq")
    
    remove_files(remove_list)
    
    

if __name__ == '__main__':
    main()
