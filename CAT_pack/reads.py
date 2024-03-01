#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 11:07:07 2020

@author: tina
"""

import os
import subprocess
import argparse
import datetime
import multiprocessing
import sys
import decimal
import gzip
from Bio import SeqIO
from io import StringIO
import json

import about
import check
import shared
import tax



def parse_arguments():
    parser = argparse.ArgumentParser(
            prog='CAT reads',
            description='Run Read Annotation Tool (RAT).',
            # @Tina: Should the explanation on how to use RAT be in the description? Right now the 'Run Read Annotation Tool' is at a little weird position. Also if you forget one argument (say -c) you'll not only get the help message but also the entire usage which is a little overkill.
            # @Tina: I do like this explanation of complete and partial workflows... Maybe we should add it to the other CAT pack as well?
            usage='CAT reads -c -t [options] [-h / --help]\n\n'
            'Complete RAT workflow (perform read mapping, run CAT, BAT, and RAT): '
            'Supply contigs, reads, database folder, taxonomy folder, and bin folder.\n\n'
            'Partial workflows:\n'
            'If you have already mapped your reads, you can supply the sorted mapping '
            'file and no read mapping will be performed.\n'
            'If you have already run CAT and/or BAT, you can supply the output '
            'files (contig2classification, bin2classification) and the path to the '
            'taxonomy folder instead.\n'
            'If you prefer not to use bin classification, do not supply the path '
            'to a bin folder.',
            add_help=False)
    
    required = parser.add_argument_group('Required arguments')
    shared.add_argument(required, 'contigs_fasta', True)
    shared.add_argument(required, 'taxonomy_folder', True)
    shared.add_argument(required, 'mode', True)
    
    optional = parser.add_argument_group('Optional arguments')
    shared.add_argument(optional, 'out_prefix', False, default='./out.RAT')
    shared.add_argument(optional, 'read_file1', False)
    shared.add_argument(optional, 'read_file2', False)
    shared.add_argument(optional, 'bam_file1', False)
    shared.add_argument(optional, 'bam_file2', False)
    shared.add_argument(optional, 'alignment_unmapped', False)
    shared.add_argument(optional, 'bin_folder', False)
    shared.add_argument(optional, 'bin_suffix', False)
    shared.add_argument(optional, 'contig2classification', False)
    shared.add_argument(optional, 'bin2classification', False)
    shared.add_argument(optional, 'read2classification', False)
    shared.add_argument(optional, 'unmapped2classification', False)

    shared.add_argument(optional, 'mapping_quality', False, default=2)
    shared.add_argument(optional, 'read_classification', False)
    shared.add_argument(optional, 'path_to_bwa', False, default='bwa')
    shared.add_argument(optional, 'path_to_samtools', False, default='samtools')
    shared.add_argument(optional, 'nproc', False,
            default=multiprocessing.cpu_count())
    shared.add_argument(optional, 'force', False)
    shared.add_argument(optional, 'quiet', False)
    shared.add_argument(optional, 'verbose', False)
    shared.add_argument(optional, 'no_log', False)
    shared.add_argument(optional, 'help', False)
    
    CAT_args = parser.add_argument_group('CAT/BAT-specific arguments')
    shared.add_argument(CAT_args, 'database_folder', False)
    shared.add_argument(CAT_args, 'r', False, default=decimal.Decimal(10))
    shared.add_argument(CAT_args, 'f', False, default=decimal.Decimal(0.5))
    shared.add_argument(CAT_args, 'path_to_prodigal', False,
            default='prodigal')
    shared.add_argument(CAT_args, 'path_to_diamond', False, default='diamond')
    shared.add_argument(CAT_args, 'no_stars', False)
    shared.add_argument(CAT_args, 'IkwId', False)
    
    dmnd_args = parser.add_argument_group('DIAMOND specific optional arguments')
    shared.add_all_diamond_arguments(dmnd_args, add_nproc=False)
                          
    (args, extra_args) = parser.parse_known_args()
    extra_args = [arg for (i, arg) in enumerate(extra_args) if
        (i, arg) != (0, 'reads')]
    if len(extra_args) > 0:
        sys.exit('error: too many arguments supplied:\n{0}'.format(
            '\n'.join(extra_args)))

    # Add extra arguments.
    shared.expand_arguments(args)
    
    if not args.contigs_fasta:
        sys.exit('error: no contigs_fasta supplied.')
    
    if not args.bam_file1 and not args.read_file1:
        sys.exit('error: you have to supply either a sorted bam file or '
                 'two paired-end read files!')
    if not args.contig2classification and not args.database_folder:
        sys.exit('error: please provide either a contig2classification file '
                 'or the path to the CAT database_folder!')
    if args.bin_folder and not args.bin2classification and not args.database_folder:
        sys.exit('error: please provide either a bin2classification file or '
                 'the path to the CAT database_folder for bin classification!')
    if args.bin2classification and not args.bin_folder:
        sys.exit('error: bin2classification file but no bin_folder supplied. '
                 'RAT requires -b/--bin_folder to include bins in annotation!')
    
    
    return args





def run():
    args = parse_arguments()

    message = '# CAT v{0}.\n'.format(about.__version__)
    shared.give_user_feedback(message, args.log_file, args.quiet,
        show_time=False)
    message = '# Running command: {0}\n\n'.format(' '.join(sys.argv[1:]))
    shared.give_user_feedback(message, args.log_file, args.quiet,
        show_time=False)
    
    # Checks
    
    
    
    errors = []
    
    errors.append(
            check.check_input_file(args.contigs_fasta, args.log_file, args.quiet))
    
    if args.read_file1:
        errors.append(
            check.check_input_file(args.read_file1, args.log_file, args.quiet))
    if args.read_file2:
        errors.append(
            check.check_input_file(args.read_file2, args.log_file, args.quiet))
    if args.bam_file1:
        errors.append(
            check.check_input_file(args.bam_file1, args.log_file, args.quiet))
    if args.contig2classification:
        errors.append(
            check.check_input_file(args.contig2classification, args.log_file, args.quiet))
    if args.bin2classification:
        errors.append(
            check.check_input_file(args.bin2classification, args.log_file, args.quiet))
    
    
    
    if not args.force:
        if args.read_file1:
            outf_bwamem='{0}.{1}.bwamem'.format(args.out_prefix+'.'+
                                        os.path.split(args.contigs_fasta)[-1], 
                                        os.path.split(args.read_file1)[-1])
            print(outf_bwamem)
            errors.append(
                    check.check_output_file(
                        outf_bwamem, args.log_file, args.quiet))
            errors.append(
                    check.check_output_file(
                        outf_bwamem + '.bam', args.log_file, args.quiet))
            errors.append(
                    check.check_output_file(
                        outf_bwamem + '.sorted', args.log_file, args.quiet))
        if not args.contig2classification:
            outf_cat_protein_faa='{0}.CAT.predicted_proteins.faa'.format(args.out_prefix)
            outf_cat_c2c='{0}.CAT.contig2classification.txt'.format(args.out_prefix)
            outf_cat_alignment='{0}.CAT.alignment.diamond'.format(args.out_prefix)
            errors.append(
                    check.check_output_file(
                        outf_cat_protein_faa, args.log_file, args.quiet))
            errors.append(
                    check.check_output_file(
                        outf_cat_c2c, args.log_file, args.quiet))
            errors.append(
                    check.check_output_file(
                        outf_cat_alignment, args.log_file, args.quiet))
        if args.bin_folder and not args.bin2classification:
            outf_bat_protein_faa='{0}.BAT.{1}.predicted_proteins.faa'.format(args.out_prefix,
                                                                         'concatenated')
            outf_bat_b2c='{0}.BAT.bin2classification.txt'.format(args.out_prefix)
            outf_bat_alignment='{0}.BAT.{1}.alignment.diamond'.format(args.out_prefix,
                                                                         'concatenated')
            errors.append(
                    check.check_output_file(
                        outf_bat_protein_faa, args.log_file, args.quiet))
            errors.append(
                    check.check_output_file(
                        outf_bat_b2c, args.log_file, args.quiet))
            errors.append(
                    check.check_output_file(
                        outf_bat_alignment, args.log_file, args.quiet))

        
    errors.append(
            check.check_samtools_binaries(
                args.path_to_samtools, args.log_file, args.quiet))
    if not args.bam_file1:
        errors.append(
            check.check_bwa_binaries(
                args.path_to_bwa, args.log_file, args.quiet))
    if not args.contig2classification or (args.bin_folder and not 
                                          args.bin2classification):
        errors.append(
            check.check_prodigal_binaries(
                args.path_to_prodigal, args.log_file, args.quiet))
        errors.append(
            check.check_diamond_binaries(
                args.path_to_diamond, args.log_file, args.quiet))
    
    if True in errors:
        sys.exit(1)
    
    if args.read_file1:
        reads_files=[args.read_file1]
        if args.read_file2:
            reads_files.append(args.read_file2)
        else:
            message = (
                'WARNING: only one read file supplied! Currently RAT does not '
                'support interlaced read files. If you are working with '
                'paired-end reads, please provide a reverse read-file!' )
            shared.give_user_feedback(message, args.log_file, args.quiet, 
                                      show_time=False)
    
    # First: run bwa mem, samtools view and samtools sort if there is no bam file
    bam_files=[]
    if not args.bam_file1:
        message = (
                '\n'
                'RAT is running. Mapping reads against assembly with bwa mem.\n')
        shared.give_user_feedback(message, args.log_file, args.quiet,
                show_time=False)
            
        

        shared.run_bwa_mem(args.path_to_bwa, args.path_to_samtools,
                              args.contigs_fasta, reads_files, args.out_prefix,
                              args.nproc, args.log_file)

            
            
        bam_files.append('{0}.{1}.bwamem.sorted'.format(args.out_prefix+'.'
                        ''+os.path.split(args.contigs_fasta)[-1], 
                        os.path.split(args.read_file1)[-1]))
    else:
        bam_files.append(args.bam_file1)
        if args.bam_file2:
            bam_files.append(args.bam_file2)
            message = (
                'WARNING: you provided two bam files. Consider mapping forward '
                'and reverse reads together or repeating the mapping step with '
                'RAT!')
            shared.give_user_feedback(message, args.log_file, args.quiet, 
                                      show_time=False)

    
    # Process bin folder
    contig2bin={}
    
    if args.bin_folder:
        message = 'Bin folder supplied. Processing bin folder.'
        shared.give_user_feedback(message, args.log_file, args.quiet,
                show_time=True)

        bins=process_bin_folder(args.bin_folder, args.bin_suffix)
        contig2bin=invert_bin_dict(bins)


    # Run BAT on folder if bin folder but not BAT_file is supplied
    # Or process BAT file if it is supplied
    b2c={}
    
    if args.bin_folder and args.bin2classification:
        message = (
                'bin2classification file supplied. Processing bin '
                'classifications.')
        shared.give_user_feedback(message, args.log_file, args.quiet,
                show_time=True)
        b2c=process_CAT_table(args.bin2classification, args.nodes_dmp, 
                              args.log_file, args.quiet)
    
    elif args.bin_folder and not args.bin2classification:
        message = ('No bin2classification file supplied. Running BAT on bin '
                   'folder.')
        shared.give_user_feedback(message, args.log_file, args.quiet,
                show_time=True)

        shared.run_BAT(args.bin_folder, args.database_folder, args.taxonomy_folder,
                       args.log_file, args.quiet, args.nproc, args.f, args.r, 
                       args.path_to_prodigal, args.path_to_diamond, args.out_prefix, 
                       args.bin_suffix)
        b2c=process_CAT_table('{0}.BAT.bin2classification.txt'.format(args.out_prefix), 
                              args.nodes_dmp, args.log_file, args.quiet)
        
    elif not args.bin_folder:
        message = 'No bin folder supplied. No bin classification will be made.'
        
        shared.give_user_feedback(message, args.log_file, args.quiet,
                show_time=False)
    
    
    # Run CAT or process CAT output file
    if args.contig2classification:
        message = (
                'contig2classification file supplied. Processing contig '
                'classifications.')
        shared.give_user_feedback(message, args.log_file, args.quiet,
                show_time=True)
        c2c=process_CAT_table(args.contig2classification, args.nodes_dmp, 
                              args.log_file, args.quiet)
    else:
        message = (
                'No contig2classification file supplied. Running CAT on '
                'contigs.')
        shared.give_user_feedback(message, args.log_file, args.quiet,
                show_time=True)
        shared.run_CAT(args.contigs_fasta, args.database_folder, 
                       args.taxonomy_folder, args.log_file, args.quiet, 
                       args.nproc, args.f, args.r, args.path_to_prodigal, 
                       args.path_to_diamond, args.out_prefix)
        c2c=process_CAT_table('{0}.CAT.contig2classification.txt'
                              ''.format(args.out_prefix), 
                              args.nodes_dmp, args.log_file, args.quiet)
    
    
    # Process BAM files and grab unmapped reads
    message = 'Processing mapping file(s).\n'
    shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=True)
    # @Tina: we could just supply the first argument as a list, that will save you the conditional.
    if len(bam_files)==2:
        reads, unmapped_reads, sum_of_reads, paired = process_bam_file(bam_files[0], 
                                                                       bam_files[1], 
                                                                       mapping_quality=args.mapping_quality)
    else:
        reads, unmapped_reads, sum_of_reads, paired = process_bam_file(bam_files[0], 
                                                                       mapping_quality=args.mapping_quality)
    
    
    # Get lengths of the contigs from contigs file    
    contig_length_dict, sum_of_nucleotides = get_contig_lengths(args.contigs_fasta)
    # make contig_dict
    contigs, reads, max_primary=make_contig_dict(reads, bam_files, paired, contig2bin)


    # If direct_mapping is chosen, map unclassified contigs and unmapped reads
    # against NR
    
    if args.mode=='direct_mapping':
        setattr(args,'read2classification',True)
        message = ('Chosen mode: {0}. Classifying unclassified contigs and'
                    ' unmapped reads with diamond.'.format(args.mode))
        shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=True)
        
        
        # unmapped_reads is a dictionary that stores the unmapped forward and
        # reverse reads and looks like this:
            # unmapped_reads={'fw': [read_id1, read_id2, read_id3], 
            #                 'rev': [read_id1, read_id2, read_id3]}
        # I need to mark all the read ids as fw or rev, pull out the sequences,
        # and then write the input fasta as one document for all reads and 
        # contigs I want to classify
        message = ('Grabbing unmapped and unclassified sequences...')
        shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=True)
        all_unclassified=list()
        
        uncl_unm_fasta='{0}.unclassified_unmapped.fasta'.format(args.out_prefix)
        unclassified_contigs = get_unclassified_contigs(contig2bin,
                                                        c2c, b2c)
        all_unclassified+=list(unclassified_contigs)
        make_unclassified_seq_fasta(args.contigs_fasta, unclassified_contigs,
                                    uncl_unm_fasta, 'fasta', 'w')
        
        message = ('Contigs written! Appending forward reads...')
        shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=True)
        make_unclassified_seq_fasta(reads_files[0], unmapped_reads['fw'],
                                    uncl_unm_fasta, 'fastq', 'a','_1')
        all_unclassified+=['{}_1'.format(i) for i in unmapped_reads['fw']]
        
        message = ('Appending reverse reads...')
        shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=True)
        make_unclassified_seq_fasta(reads_files[1], unmapped_reads['rev'],
                                    uncl_unm_fasta, 'fastq', 'a','_2')
        all_unclassified+=['{}_2'.format(i) for i in unmapped_reads['rev']]
        # Run diamond on unclassified contigs and unmapped reads
        if not args.alignment_unmapped and not args.unmapped2classification: 

            setattr(args,
                    'alignment_file',
                    '{0}.unclassified_unmapped.alignment.diamond'.format(args.out_prefix))
            shared.run_diamond(args, blast='blastx', prot_fasta=uncl_unm_fasta,
                                top=11)
        else:
            setattr(args,
                    'alignment_file',
                    args.alignment_unmapped)
        # Now, the diamond alignment has to be parsed with a CAT function
        if not args.unmapped2classification:
    
            taxid2parent, taxid2rank = tax.import_nodes(args.nodes_dmp, 
                                                        args.log_file,
                                                        args.quiet)
            shared.explore_database_folder(args)
            seq2hits, all_hits = shared.parse_tabular_alignment(args.alignment_file,
                                                                decimal.Decimal(1-args.r),
                                                                args.log_file,
                                                           args.quiet)
            fastaid2LCAtaxid = tax.import_fastaid2LCAtaxid(args.fastaid2LCAtaxid_file, 
                                                           all_hits, args.log_file,
                                                           args.quiet)
            taxids_with_multiple_offspring = tax.import_taxids_with_multiple_offspring(
                args.taxids_with_multiple_offspring_file,
                args.log_file,
                args.quiet)
            message = ('Finding lineages for unclassified/unmapped sequences...')
            shared.give_user_feedback(message, args.log_file, args.quiet,
                show_time=True)
            write_unmapped2classification(seq2hits, 
                                          all_unclassified,
                                          fastaid2LCAtaxid, 
                                          args.out_prefix, 
                                          taxid2parent,
                                          taxid2rank,
                                          taxids_with_multiple_offspring,
                                          args.no_stars)
            u2c=process_CAT_table('{0}.unmapped2classification.txt'.format(args.out_prefix), 
                                  args.nodes_dmp, 
                                  args.log_file, 
                                  args.quiet)
        else:
            u2c=process_CAT_table(args.unmapped2classification, 
                                  args.nodes_dmp, 
                                  args.log_file, 
                                  args.quiet)
    else:
        u2c=None
            
    message = 'Writing output tables.'
    shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=True)    
    make_tax_table(c2c,
                   contigs, 
                   unmapped_reads,
                   contig_length_dict,
                   contig2bin, 
                   b2c, 
                   sum_of_nucleotides, 
                   sum_of_reads, 
                   args.out_prefix,
                   u2c)
    
    make_bin_table(contig2bin, 
                   contigs, 
                   contig_length_dict, 
                   sum_of_reads, 
                   args.out_prefix)
    
       
    
    # Finally: Optional: Classify reads and write per_read table
    if args.read2classification:
        reads=classify_reads(reads, c2c, b2c, u2c)
        write_read_table(reads, args.out_prefix, max_primary)
        # print("Writing read dictionary to file...")
        # json.dump(reads, open('/home/tina/Documents/RAT/debug_direct_mapping/'
        #                       '20211104_read_dict.json', 'w'))
        # print("Written!")
    
    message = (
            '\n-----------------\n\n'
            '[{0}] RAT is done!'.format(
                datetime.datetime.now()))
    shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=False)
    return





def write_unmapped2classification(seq2hits,
                                  all_seqs,
                                  fastaid2LCAtaxid, 
                                  out_prefix,
                                  taxid2parent, 
                                  taxid2rank, 
                                  taxids_with_multiple_offspring,
                                  no_stars):
    
    
    
    with open('{0}.unmapped2classification.txt'.format(out_prefix), 'w') as outf:
        n_classified_contigs = 0
        outf.write(
                '# sequence\tclassification\treason\tlineage\tlineage scores\n')
        
        for seq in sorted(all_seqs):
            if seq not in seq2hits:
                outf.write('{0}\tno taxid assigned\tno hits found\n'.format(
                    seq))
                continue
            
            n_hits = len(seq2hits[seq])
            LCAs_ORFs = []
            
            
            (taxid,
                    top_bitscore) = tax.find_LCA_for_ORF(
                            seq2hits[seq], fastaid2LCAtaxid, taxid2parent)
             
            if not taxid.startswith('no taxid found'):
                lineage = tax.find_lineage(taxid, taxid2parent)
    
                if not no_stars:
                    lineage = tax.star_lineage(
                            lineage, taxids_with_multiple_offspring)
                               
            LCAs_ORFs.append((taxid, top_bitscore))
                    
            if len(LCAs_ORFs) == 0:
                outf.write('{0}\tno taxid assigned\t'
                        'no hits to database\n'.format(seq))
    
                continue
            
            (lineages,
                    lineages_scores,
                    based_on_n_ORFs) = tax.find_weighted_LCA(
                            LCAs_ORFs, taxid2parent, 0.5)
             
            if lineages == 'no ORFs with taxids found.':
                outf.write('{0}\tno taxid assigned\t'
                        'hits not found in taxonomy files\n'.format(seq))
    
                continue
            
            if lineages == 'no lineage whitelisted.':
                outf.write(
                        '{0}\tno taxid assigned\t'
                        'no lineage reached minimum bit-score support\n'
                        ''.format(seq))
    
                continue
        
            n_classified_contigs += 1

            for (i, lineage) in enumerate(lineages):
                if not no_stars:
                    lineage = tax.star_lineage(
                            lineage, taxids_with_multiple_offspring)
                scores = ['{0:.2f}'.format(score) for
                        score in lineages_scores[i]]
                
                if len(lineages) == 1:
                    # There is only one classification.
                    outf.write(
                            '{0}\t'
                            'taxid assigned\t'
                            'based on {1} hits\t'
                            '{2}\t'
                            '{3}\n'.format(
                                seq,
                                n_hits,
                                ';'.join(lineage[::-1]),
                                ';'.join(scores[::-1])))
     
    return




def classify_reads(read_dict, contig2classification, bin2classification, unmapped2classification):
    worked=0
    for read in read_dict:
        # If the read is mapped to a contig:
        if read_dict[read]['contig']!=[]:
            for r in range(len(read_dict[read]['contig'])):
                # If there is bin classification, store it
                if bin2classification:
                    try:
                        if read_dict[read]['bin'][r]=='unbinned':
                            read_dict[read]['taxon_bin'].append('')
                        else:
                            read_dict[read]['taxon_bin'].append(';'.join(bin2classification[read_dict[read]['bin'][r]][0]))
                        worked+=1
                    except IndexError:
                        print(r)
                        sys.exit('Worked: {}'.format(worked))
                
                # If there is a contig classification, store it
                contig=read_dict[read]['contig'][r]
                if (contig2classification[contig]==[] or 
                    len(contig2classification[contig][0])<2 or 
                    (len(contig2classification[contig][0])==2 and 
                    contig2classification[contig][0][1]=='131567')):
                    read_dict[read]['taxon_contig'].append('')
                else:
                    read_dict[read]['taxon_contig'].append(';'.join(contig2classification[contig][0]))
                    
                # If there is no contig classification, but a dm classification, store it
                if contig in unmapped2classification:
                    if (unmapped2classification[contig]==[] or 
                    len(unmapped2classification[contig][0])<2 or 
                    (len(unmapped2classification[contig][0])==2 and 
                    unmapped2classification[contig][0][1]=='131567')):
                        read_dict[read]['taxon_contig_dm'].append('')
                    else:
                        read_dict[read]['taxon_contig_dm'].append(';'.join(unmapped2classification[contig][0]))
        # If the read is not mapped:
        else:
            fw='{}_1'.format(read)
            rev='{}_2'.format(read)
            if fw in unmapped2classification:
                if (unmapped2classification[fw]==[] or
                len(unmapped2classification[fw][0])<2 or 
                (len(unmapped2classification[fw][0])==2 and 
                 unmapped2classification[fw][0][1]=='131567')):
                    read_dict[read]['taxon_read_dm'].append('')
                else:
                    read_dict[read]['taxon_read_dm'].append('fw: '+';'.join(unmapped2classification[fw][0]))
            if rev in unmapped2classification:
                if (unmapped2classification[rev]==[] or
                len(unmapped2classification[rev][0])<2 or 
                (len(unmapped2classification[rev][0])==2 and 
                 unmapped2classification[rev][0][1]=='131567')):
                    read_dict[read]['taxon_read_dm'].append('')
                else:
                    read_dict[read]['taxon_read_dm'].append('rev: '+';'.join(unmapped2classification[rev][0]))
    
    return read_dict



def write_read_table(read_dict, sample_name, max_primary):
    # Takes the read dictionary with bin lineage and contig lineage and writes
    # annotation for each read and each tier of annotation into a file
    read_table=sample_name + '.read2classification.txt'
    
    with open(read_table, 'w') as outf:
        outf.write('## command: {}\n'.format(' '.join(sys.argv)))
        outf.write('# read\tclassification\tbin classification\t'
                   'contig classification\tdirect_mapping_contig\t'
                   'direct_mapping_read\n'
                   #'[aligner] unclassified contig\t[aligner] unclassified read\n'
                   )
        for read in read_dict:

            for c in range(len(read_dict[read]['contig'])):
                try:
                    b_taxon=read_dict[read]['taxon_bin'][c]
                except IndexError:
                    b_taxon=''
                try:    
                    c_taxon=read_dict[read]['taxon_contig'][c]
                except IndexError:    
                    c_taxon=''
                if b_taxon=='' and c_taxon=='':
                    try:
                        c_dm_taxon=read_dict[read]['taxon_contig_dm'][c]
                        r_dm_taxon=''
                    except IndexError:
                        c_dm_taxon=''
                        r_dm_taxon=''
                    
                else:
                    c_dm_taxon=''
                    r_dm_taxon=''
                assigned='taxid assigned'
                if b_taxon=='' and c_taxon=='' and c_dm_taxon=='':
                    assigned='no taxid assigned'
                outf.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(read,
                           assigned, b_taxon, c_taxon, c_dm_taxon,
                           r_dm_taxon))
            if len(read_dict[read]['contig'])<max_primary:
                if len(read_dict[read]['taxon_read_dm'])>0:
                    for c in range(len(read_dict[read]['taxon_read_dm'])):
                
                        try:
                            r_dm_taxon=read_dict[read]['taxon_read_dm'][c]
                            
                        except IndexError:
                            r_dm_taxon=''
                        assigned='taxid assigned'
                        if r_dm_taxon=='':
                            assigned='no taxid assigned'
                        outf.write('{0}\t{1}\t\t\t\t{2}\n'.format(read,
                           assigned, r_dm_taxon))
                else:
                    n=max_primary-len(read_dict[read]['taxon_read_dm'])
                    outf.write(n*'{0}\t{1}\n'.format(read, 'no taxid assigned'))
    return



def make_bin_table(contig2bin, 
                   contig_dict, 
                   contig_length_dict, 
                   total_reads, 
                   out_prefix):
    # Include lineage?
    bin_data={}
    bin_file='{0}.bin.reads.txt'.format(out_prefix)
    
    with open(bin_file, 'w') as op:
        op.write('## command: {}\n'.format(' '.join(sys.argv)))
        op.write('# bin\tnumber of reads\tfraction of reads\tbin length\t'
                 'corrected fraction\n')
        for contig in contig_dict:
            if contig in contig2bin:
                bin_id=contig2bin[contig]
            else:
                bin_id='unbinned'
            if not bin_id in bin_data:
                bin_data[bin_id]={
                    'n_reads': 0,
                    'n_nucleotides': 0
                    }
            # if contig is in contig_dict, then reads have mapped to it, otherwise
            # no reads have mapped to it
            try:
                bin_data[bin_id]['n_reads']+=contig_dict[contig]
                bin_data[bin_id]['n_nucleotides']+=contig_length_dict[contig]
            except KeyError:
                pass
            
        divisor=get_proper_fraction(bin_data, total_reads)
        for bin_id in bin_data:
            read_fraction=bin_data[bin_id]['n_reads']/total_reads
            if bin_data[bin_id]['n_nucleotides']!=0:
                corrected_frac=(read_fraction/bin_data[bin_id]['n_nucleotides'])/divisor
            else:
                corrected_frac=0
            op.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(bin_id, 
                                                bin_data[bin_id]['n_reads'], 
                                                read_fraction, 
                                                bin_data[bin_id]['n_nucleotides'], 
                                                corrected_frac))
    return





def make_tax_table(c2c_dict, 
                   contig_dict,
                   unmapped_reads,
                   contig_length_dict, 
                   contig2bin, 
                   bin2classification, 
                   sum_of_nucleotides, 
                   total_reads, 
                   sample_name,
                   u2c_dict=None):
    """
    Makes a dictionary that keeps all the necessary info in the same place
    Writes contig table
    Writes table with all classifications including subspecies
    """
    RAT_contig_header=('# contig\tnumber of reads\tfraction of reads\t'
                       'contig_length\tAverage Coveraget\tlineage\t'
                       'lineage ranks\tlineage scores\n')
                       
    RAT_file_header=('# lineage\tnumber of reads\tfraction of reads\t'
                     'taxon length\tlineage ranks\n')
    
    taxon_dict={}
    unclassified={
            'contigs': set(),
            'n_reads': 0,
            'n_nucleotides': 0
            }

    for contig in contig_dict:
        
        # If the contig name is '*', that means that a read is unmapped. 
        # If direct_mapping is enabled, unmapped reads have been run through
        # diamond and an unmapped2classification dictionary has been supplied
        # to this function. The reads get assigned a taxon if diamond spit out
        # a taxon, they get counted as unmapped otherwise.
        
        if contig=='*':
            if u2c_dict:
                unmapped=0
                tmp=['fw', 'rev']
                for direction in tmp:
                    for read in unmapped_reads[direction]:
                        u2c_read_id='{0}_{1}'.format(read, 
                                                     tmp.index(direction)+1)
                        if u2c_read_id in u2c_dict:
                            if (u2c_dict[u2c_read_id]==[] or len(u2c_dict[u2c_read_id][0])<2 
                                or (len(u2c_dict[u2c_read_id][0])==2 and 
                                u2c_dict[u2c_read_id][0][1]=='131567')):
                                unmapped+=1
                            else:
                                taxon=';'.join(u2c_dict[u2c_read_id][0])
                                if taxon not in taxon_dict:
                                    taxon_dict[taxon]={
                                            'n_reads': 0,
                                            'n_nucleotides': 0,
                                            'ranks': u2c_dict[u2c_read_id][1]
                                            }
                                taxon_dict[taxon]['n_reads']+=1
                                # taxon_dict[taxon]['n_nucleotides']+=contig_length_dict[contig]
            else:
                unmapped=contig_dict[contig]
        
        # If the contig is in contig2bin, that means that it is binned and therefore
        # will be assigned the lineage of the bin and not the contig
        elif contig in contig2bin:
            taxon=';'.join(bin2classification[contig2bin[contig]][0])
            if taxon not in taxon_dict:
                taxon_dict[taxon]={
                        'n_reads': 0,
                        'n_nucleotides': 0,
                        'ranks': bin2classification[contig2bin[contig]][1]
                        }
            # add number of reads that map to the contig to the taxon and
            # add nucleotides of the contig to the taxon
            taxon_dict[taxon]['n_reads']+=contig_dict[contig]
            taxon_dict[taxon]['n_nucleotides']+=contig_length_dict[contig]
            
        # If C2C[contig] is an empty list, it's unclassified
        # If the length of the list of classified tax_ids is shorter than 2, the
        # contig is not classified at superkingdom level and therefore, for our
        # purpose unclassified
        # We check u2c to see whether it was classified in the direct_mapping
        # step, assign the lineage from u2c or unclassified depending on
        # result
        elif (c2c_dict[contig]==[] or len(c2c_dict[contig][0])<2 
              or (len(c2c_dict[contig][0])==2 and 
              c2c_dict[contig][0][1]=='131567')):
            if u2c_dict:
                if contig not in u2c_dict or (u2c_dict[contig]==[] or len(u2c_dict[contig][0])<2 
                  or (len(u2c_dict[contig][0])==2 and 
                  u2c_dict[contig][0][1]=='131567')):
                    unclassified['n_reads']+=contig_dict[contig]
                    unclassified['n_nucleotides']+=contig_length_dict[contig]
                else:
                    taxon=';'.join(u2c_dict[contig][0])
                    if taxon not in taxon_dict:
                        taxon_dict[taxon]={
                                'n_reads': 0,
                                'n_nucleotides': 0,
                                'ranks': u2c_dict[contig][1]
                                }
                    taxon_dict[taxon]['n_reads']+=contig_dict[contig]
                    taxon_dict[taxon]['n_nucleotides']+=contig_length_dict[contig]
                
            
        
        # The other possibility is that the contig is unbinned, but classified
        else:
            taxon=';'.join(c2c_dict[contig][0])
            if taxon not in taxon_dict:
                taxon_dict[taxon]={
                        'n_reads': 0,
                        'n_nucleotides': 0,
                        'ranks': c2c_dict[contig][1]
                        }
            
            # add number of reads that map to the contig to the taxon and
            # add nucleotides of the contig to the taxon
            taxon_dict[taxon]['n_reads']+=contig_dict[contig]
            taxon_dict[taxon]['n_nucleotides']+=contig_length_dict[contig]
    
    unmapped_line='{0}\t{1}\t{2}\n'.format('unmapped', unmapped,
                                           unmapped/total_reads)
    unclassified_line='{0}\t{1}\t{2}\t{3}\t'.format('unclassified',
                       unclassified['n_reads'],
                       unclassified['n_reads']/total_reads,
                       unclassified['n_nucleotides'])
    
    
    # write the output tables for contigs:
    with open('{0}.contig.abundance.txt'.format(sample_name), 'w') as outf_contig:
        outf_contig.write('## command: {}\n'.format(' '.join(sys.argv)))
        outf_contig.write(RAT_contig_header)
        # @Tina: not to work with decimals instead of floats.
        # get the divisor to make corrected fraction add up to 1

        divisor=get_proper_fraction(contig_dict, 
                                    total_reads, 
                                    unclassified, 
                                    contig_length_dict)
        outf_contig.write(unmapped_line)
        for contig in contig_dict:
            if not contig=='*':
                read_frac=contig_dict[contig]/total_reads
                corrected_frac=(read_frac/contig_length_dict[contig])/divisor
                # If contig is binned, then get the lineage of the bin and not the contig
                try:
                    if contig in contig2bin:
                        lineage=';'.join(bin2classification[contig2bin[contig]][0])
                        rank_list=';'.join(bin2classification[contig2bin[contig]][1])
                        scores=';'.join(bin2classification[contig2bin[contig]][2])
                        
                    else:
                        lineage=';'.join(c2c_dict[contig][0])
                        rank_list=';'.join(c2c_dict[contig][1])
                        scores=';'.join(c2c_dict[contig][2])
                        
                except IndexError:
                    lineage=''
                    rank_list=''
                    scores=''
                outf_contig.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(
                        contig, contig_dict[contig], read_frac, 
                        contig_length_dict[contig], corrected_frac, lineage, 
                        rank_list, scores))
    
    
    # write output table for complete lineage
    with open('{0}.complete.abundance.txt'.format(sample_name), 'w+') as outf_complete:
        outf_complete.write('## command: {}\n'.format(' '.join(sys.argv)))
        outf_complete.write(RAT_file_header)
        outf_complete.write(unmapped_line)
        outf_complete.write(unclassified_line+'\n')

        # if unclassified['n_reads']!=0:
        #     outf_complete.write('{0}\n'.format(((unclassified['n_reads']/total_reads)/unclassified['n_nucleotides'])/divisor))
        # else:
        #     outf_complete.write('0\n')
        
        for taxon in taxon_dict:
            read_frac=taxon_dict[taxon]['n_reads']/total_reads

            outf_complete.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(
                    taxon, 
                    taxon_dict[taxon]['n_reads'],
                    read_frac, 
                    taxon_dict[taxon]['n_nucleotides'],
                    # add average genome size
                    ';'.join(taxon_dict[taxon]['ranks'])))
    
    
    return 




def get_proper_fraction(rank_dict, 
                        total_reads, 
                        unclassified={}, 
                        contig_length_dict={}):
    divisor=0
    if unclassified and unclassified['n_reads']!=0:
        divisor+=(unclassified['n_reads']/total_reads)/unclassified['n_nucleotides']
    try:
        for taxon in rank_dict:
            if not taxon=='*':
                if rank_dict[taxon]['n_nucleotides']!=0:
                    divisor+=(rank_dict[taxon]['n_reads']/total_reads)/rank_dict[taxon]['n_nucleotides']
    except TypeError:
        for contig in rank_dict:
            if not contig=='*':
                divisor+=(rank_dict[contig]/total_reads)/contig_length_dict[contig]
    return divisor



def process_CAT_table(CAT_output_table, 
                      nodes_dmp, 
                      log_file, 
                      quiet):
    """
    takes CAT contig2classification or bin2classification table and taxonomy 
    nodes file as input; returns a dictionary that contains the complete 
    lineage for each contig/bin:
    (list of tuples of the format [(tax_id1,tax_id2...), 
                                   (rank1,rank2,...), 
                                   (score1,score2...)])
    bin id contains the file extension!
    """
    c2c={}
    
    # get rank for each taxid
    taxid2rank=tax.import_nodes(nodes_dmp, log_file, quiet)[1]
    
    with open(CAT_output_table, 'r') as contigs:
        for contig in contigs:
            
            # if it is not the first line, get contig id, lineage and lineage scores
            if not contig.startswith('# '):
                contig_id=contig.strip().split('\t')[0]
                c2c[contig_id]=[]
                
                # Build a list of tuples for the lineage:
                # (tax_id1,tax_id2...), (rank1,rank2,...), (score1,score2...)]
                if len(contig.split('\t'))>3:
                    tax_ids=contig.split('\t')[3].strip().split(';')
                    scores=contig.split('\t')[4].strip().split(';')
                    c2c[contig_id]=[tuple([t for t in tax_ids]),
                                   tuple([taxid2rank[t.strip('*')] for t in tax_ids]),
                                   tuple([s for s in scores])]
                                        
    return c2c
            




def process_bam_file(BAM_fw_file, BAM_rev_file=False, path_to_samtools='samtools', mapping_quality=2):
    """
    # function that processes the bam file(s) by picking out the primary alignments
    # and storing the contigs that the forward and reverse reads map to. If only
    # one file is given, checks whether it's paired or not before it counts the
    # unmapped reads.
    # Keeps track of unmapped reads: If paired read file, then if checks the flag
    # for whether the unmapped read is the first or the second in the pair; if
    # it's single reads, it just fills up the fw key of the dictionary, and if
    # there is a reverse file, then the reverse file only adds to the rev key.
    """
    read_dict = {}
    unmapped_reads = {'fw': set(), 'rev': set()}
    
    # Open alignment file with samtools to avoid dependency on pysam
    cmd = [path_to_samtools, 'view', BAM_fw_file]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    # Check whether it is a paired-end BAM file by examining the flags
    paired = False

    # Build read dictionary with list of 0-2 contigs that each read/mate maps to
    for mapping in proc.stdout:
        mapping = mapping.decode("utf-8").rstrip().split('\t')
        read_id, flag, contig, score = mapping[0], int(mapping[1]), mapping[2], int(mapping[4])
        
        if not paired and bin(flag)[-1] == '1':
            paired = True
        
        if read_id not in read_dict:
            read_dict[read_id] = {
                'contig': [],
                'bin': [],
                'taxon_bin': [],
                'taxon_contig': [],
                'taxon_contig_dm': [],
                'taxon_read_dm': [],
            }
        
        # Check if the score is above the threshold and the alignment is not supplementary
        if score >= mapping_quality and len(bin(flag)) < 14:
            read_dict[read_id]['contig'].append(contig)
        
        # Store unmapped reads as forward and reverse
        elif contig == '*' or (score < mapping_quality and len(bin(flag)) < 14):
            # Check if read is the "first in pair"
            if paired:
                if len(bin(flag)) > 8 and bin(flag)[-7] == '1':
                    unmapped_reads['fw'].add(read_id)
                else:
                    unmapped_reads['rev'].add(read_id)
            else:
                unmapped_reads['fw'].add(read_id)
    
    
    # If a reverse alignment file is given, add the contig that the mate maps to
    if BAM_rev_file:
        cmd = [path_to_samtools, 'view', BAM_rev_file]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        for read in proc.stdout:
            read = read.decode("utf-8").rstrip().split('\t')
            read_id, flag, contig, score = read[0], int(read[1]), read[2], int(read[4])
            if score >= mapping_quality and len(bin(flag)) < 14:
                read_dict[read_id]['contig'].append(contig)
            elif contig == '*' or (score < mapping_quality and len(bin(flag)) < 14):
                unmapped_reads['rev'].add(read_id)
    
    # Calculate the sum of reads as the number of reads in the dictionary (including unmapped reads)
    sum_of_reads = 2 * len(read_dict) if paired or BAM_rev_file else len(read_dict)

    return read_dict, unmapped_reads, sum_of_reads, paired



def make_contig_dict(read_dict, 
                     bam_files, 
                     paired, 
                     contig2bin):
    """
    Takes the read dictionary that stores which read maps where and returns
    a dictionary that stores how many reads bind to each contig. It also adds
    bin annotation to the read dict.
    """
    contig_dict = {}
    contig_dict['*'] = 0
    
    # Check the maximum number of primary alignments allowed per read
    max_primary_alignment = 2 if len(bam_files) > 1 or paired else 1
    
    for read in read_dict:
        # Check if a read (pair) has too many primary alignments
        if len(read_dict[read]['contig']) > max_primary_alignment:
            message = ('Something is wrong with the alignment. Too many '
                       'mappings detected. Please check your alignment file or' 
                       ' map your reads using RAT.')
            shared.give_user_feedback(message)
            sys.exit('Mapping file is wrong format.')
        
        # Otherwise, calculate the number of unmapped reads and update the contig dictionary
        else:
            contig_dict['*'] += max_primary_alignment - len(read_dict[read]['contig'])
            for contig in read_dict[read]['contig']:
                if contig not in contig_dict:
                    contig_dict[contig] = 0
                
                contig_dict[contig] += 1
                if contig2bin:
                    if contig in contig2bin:
                        read_dict[read]['bin'].append(contig2bin[contig])
                    else:
                        read_dict[read]['bin'].append('unbinned')
                    
    return contig_dict, read_dict, max_primary_alignment

    


def make_unclassified_seq_fasta(seq_fasta, 
                                unclassified_seq_ids, 
                                unclassified_seq_fasta, 
                                f_format, 
                                f_mode, 
                                suffix=''):
    """
    Extracts the unclassified sequences from a fasta file based on the provided
    IDs and creates a new fasta file with those sequences.
    Supports gzip-compressed fasta files as well.
    The resulting fasta file contains the unclassified sequences with an 
    optional suffix in the sequence ID.
    """
    
    # If fasta is compressed fasta file
    if seq_fasta.endswith('.gz'):
        fasta_dict = {}
        fasta_file = gzip.open(seq_fasta, 'rt')
        record = ''  # Variable to store the current fasta record

        # Loop over the fasta_file and pick out the unclassified ones.
        for line in fasta_file:
            tmp = str(line)
            
            if f_format=='fastq':
                new_seq='@'
            elif f_format=='fasta':
                new_seq='>'
                
            if tmp.startswith(new_seq) and record:
                seq_id = record[1:].split()[0]

                if seq_id in unclassified_seq_ids:
                    tmp2 = SeqIO.to_dict(SeqIO.parse(StringIO(record), f_format))
                    fasta_dict.update(tmp2) 
                record = tmp
            else:
                record += tmp

        tmp2 = SeqIO.to_dict(SeqIO.parse(StringIO(record), f_format))
        seq_id = record[1:].split()[0]
        if seq_id in unclassified_seq_ids:
            fasta_dict.update(tmp2)  # Add the unclassified sequence record to the dictionary
    
    # If fasta isn't compressed
    else:
        fasta_file = open(seq_fasta)
         # Convert fasta records to dictionary
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, f_format)) 

    with open(unclassified_seq_fasta, f_mode) as outf:
        suffices = ['/1', '/2', '_1', '_2']
        for seq in unclassified_seq_ids:
            if seq in fasta_dict:
                seq = seq
            else:
                for s in suffices:
                    if seq + s in fasta_dict:
                        seq += s
            outf.write('>{0}{1}\n{2}\n'.format(fasta_dict[seq].id, suffix, fasta_dict[seq].seq))

    return



def get_unclassified_contigs(contig2bin, c2c, b2c):
    """
    If the contig is in a bin that is unclassified, or if the contig is not in
    a bin AND unclassified: put it in a set.
    """
    unclassified_contigs = []

    for c in c2c:
        classified = False

        # Check if contig is in a bin and if the bin is classified
        if c in contig2bin:
            mag_lineage = b2c.get(contig2bin[c], [])
            if len(mag_lineage) > 2 or (len(mag_lineage) == 2 and mag_lineage[1] != '131567'):
                classified = True

        # If contig is not classified in any bin, check if the contig itself is classified
        if not classified:
            lineage = c2c.get(c, [])
            if len(lineage) > 2 or (len(lineage) == 2 and lineage[1] != '131567'):
                classified = True

        # If contig is not classified in any bin or by itself, add it to the list of unclassified contigs
        if not classified:
            unclassified_contigs.append(c)

    return unclassified_contigs
        




def get_contig_lengths(contig_file):
    """
    gets the lengths of all contigs from the contig file. Works for all types
    of contig IDs, as long as it is a proper file with unique headers, because
    it checks the length of the strings.
    """
    contig_length_dict = {}
    sum_of_nucleotides = 0

    with open(contig_file, 'r') as inp:
        for line in inp:
            line = line.strip()
            if line.startswith('>'):
                contig_id = line[1:]
                contig_length_dict[contig_id] = 0
            else:
                contig_length_dict[contig_id] += len(line)
                sum_of_nucleotides += len(line)

    return contig_length_dict, sum_of_nucleotides






def process_bin_folder(path_to_folder, bin_suffix):
    """
    takes the path to a folder containing metagenomic bins and picks out all
    nucleotide fasta files. Returns a dictionary that contains the bin ids as
    keys and a set of all contig ids that belong in that bin as values. 
    Important: Only contains the binned sequences. Every contig that is not
    in the output dictionary is unbinned.
    """
    bin_dict = {}

    # Get all files in the folder
    files = os.listdir(path_to_folder)

    for file_name in files:
        if file_name.endswith(bin_suffix) and not file_name.startswith('unbinned'):
            bin_id = file_name.strip()
            bin_dict[bin_id] = set()

            # Open fasta file and split it into contigs
            with open(os.path.join(path_to_folder, file_name), 'r') as fasta_file:
                contigs = fasta_file.read().split('>')[1:]

                # Add contig ids to the set
                for contig in contigs:
                    contig_id = contig.split()[0]
                    bin_dict[bin_id].add(contig_id)

    return bin_dict




def invert_dict(bin_dict):
    
    contig2bin = {bin_value: key for key, bin_value in bin_dict.items()}

    return contig2bin



# if __name__ == '__main__':
#     sys.exit('Run \'CAT\' to run CAT, BAT, or RAT.')
