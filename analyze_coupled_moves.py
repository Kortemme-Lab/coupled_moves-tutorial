import os
import sys
import weblogolib
import corebio
import subprocess
import numpy as np
import random
import shutil
import cPickle as pickle
from klab.Reporter import Reporter
from cogent import LoadSeqs, PROTEIN, DNA, RNA
from cogent.core.alignment import DenseAlignment
from cogent.evolve.coevolution import validate_alignment, coevolve_pair_functions, coevolve_pair, sca_input_validation, sca_pair, ancestral_state_pair, ancestral_states_input_validation, validate_position, coevolve_alignment, coevolve_alignment_functions
# import multiprocessing

class Seqs:
    def __init__(self):
        self.alphabet = corebio.seq.unambiguous_protein_alphabet
        self.alphabet_as_list_shuffled = [x for x in self.alphabet]
        random.shuffle( self.alphabet_as_list_shuffled )
        self.raw_seqs = corebio.seq.SeqList( alphabet = self.alphabet )
        self.seq_length = None
        self.design_positions = None
        self.starting_seq = None
        self.name = None
        self.seq_dicts = []
        self.dir_path = None
        self.mi = {}

    def make_wt(self):
        wt_seqs = Seqs()
        wt_seqs.alphabet = self.alphabet
        wt_seqs.alphabet_as_list_shuffled = self.alphabet_as_list_shuffled
        wt_seqs.raw_seqs = corebio.seq.SeqList( alphabet = self.alphabet )
        wt_seqs.seq_length = self.seq_length
        wt_seqs.design_positions = self.design_positions
        wt_seqs.starting_seq = self.starting_seq
        wt_seqs.name = self.name + '_wt'
        wt_seqs.seq_dicts = []
        wt_seqs.dir_path = self.dir_path
        wt_seqs.raw_seqs.append(self.starting_seq)
        return wt_seqs

    @staticmethod
    def self_minus_other(s, o):
        delta_seqs = Seqs()
        delta_seqs.load_seqs_from_freqs_matrix( s.delta_other(o) )
        assert( s.design_positions == o.design_positions )
        delta_seqs.design_positions = s.design_positions
        assert( o.starting_seq == s.starting_seq )
        delta_seqs.starting_seq = s.starting_seq
        if s.name and o.name:
            delta_seqs.name = '%s-%s' % (s.name, o.name)
        return delta_seqs

    def mutations_enriched_over_other(self, other):
        delta_matrix = self.delta_other(other)
        self_matrix = self.as_matrix()
        other_matrix = other.as_matrix()
        mutations = []
        for i, delta_matrix_position in enumerate(delta_matrix):
            design_position = self.design_positions[i]
            for j, enrichment_factor in enumerate(delta_matrix_position):
                mutations.append( (enrichment_factor, design_position, self.alphabet[j], self.starting_seq[i], self_matrix[i][j], other_matrix[i][j], i) )
        mutations.sort(reverse = True)
        return mutations

    def load_seqs_from_freqs_matrix(self, matrix):
        decoy_seqs_to_make = 10000
        assert( len(matrix.shape) == 2 )
        self.seq_length = matrix.shape[0]
        assert( len(self.alphabet) == matrix.shape[1] )
        all_position_chars = []
        for pos_index in xrange(self.seq_length):
            position_chars = []
            # # Rescale to 1.0
            # if np.sum(matrix[pos_index]) > 0.0:
            #     scale_factor = np.float64(1.0) / np.sum(matrix[pos_index])
            #     rescaled_row = matrix[pos_index] * scale_factor
            # else:
            #     rescaled_row = matrix[pos_index]
            rescaled_row = matrix[pos_index]
            for char_index, char in enumerate(self.alphabet):
                if rescaled_row[char_index] > 0.0:
                    for char_num in xrange( int( np.around(np.float64(decoy_seqs_to_make) * rescaled_row[char_index]) ) ):
                        position_chars.append(char)
            random.shuffle(position_chars)
            while len(position_chars) > decoy_seqs_to_make:
                position_chars.pop()
            # random.shuffle(self.alphabet_as_list_shuffled)
            # alphabet_index = 0
            # while len(position_chars) < decoy_seqs_to_make:
            #     position_chars.append( self.alphabet_as_list_shuffled[alphabet_index] )
            #     alphabet_index = (alphabet_index + 1) % len(self.alphabet_as_list_shuffled)
            if len(position_chars) == 0:
                while len(position_chars) < decoy_seqs_to_make:
                    position_chars.append( '' )
            while len(position_chars) < decoy_seqs_to_make:
                position_chars.append( random.choice(position_chars) )
            random.shuffle(position_chars)
            all_position_chars.append(position_chars)

        for i in xrange(decoy_seqs_to_make):
            decoy_seq = ''
            for pos_index in xrange(self.seq_length):
                decoy_seq += all_position_chars[pos_index][i]
            self.raw_seqs.append(decoy_seq)

    def load_seqs_from_stats(self, stats_path):
        new_seqs = []
        seq_count = 0
        with open(stats_path, 'r') as f:
            for line in f:
                data = line.strip().split()
                seq_dict = { 'scores' : {'seq_count' : seq_count} }
                seq_count += 1
                for i in xrange(len(data)):
                    if i % 2 != 0 or i + 1 >= len(data):
                        # continue for odd numbers and if no pair
                        pass
                    elif i == 0:
                        seq_dict['scores']['total_score'] = float(data[i+1])
                    else:
                        label = data[i].strip(':')
                        value = data[i+1]
                        try:
                            float_value = float(value)
                            seq_dict['scores'][label] = float_value
                        except ValueError:
                            seq_dict[label] = value
                            if label == 'sequence':
                                if self.seq_length == None:
                                    self.seq_length = len(value)
                                else:
                                    assert( len(value) == self.seq_length )
                new_seqs.append(seq_dict)
        self.seq_dicts.extend( new_seqs )
        self.raw_seqs.extend( [d['sequence'] for d in new_seqs] )

    def filter_by_score_percentile(self, score_types, percentile_cutoff, reverse = False):
        raise Exception( 'Check if reverse should be False before using' )
        score_list = [(sum(seq_dict['scores'][score_type] for score_type in score_types), seq_dict['scores']['seq_count'], seq_dict) for seq_dict in self.seq_dicts]
        score_list.sort( reverse = reverse )
        filtered_scores = score_list[ : int( float(len(score_list)) * float(percentile_cutoff) )]
        filtered_seqs = Seqs()
        filtered_seqs.alphabet = self.alphabet
        filtered_seqs.alphabet_as_list_shuffled = self.alphabet_as_list_shuffled
        filtered_seqs.raw_seqs = corebio.seq.SeqList( alphabet = self.alphabet )
        filtered_seqs.seq_length = self.seq_length
        filtered_seqs.design_positions = self.design_positions
        filtered_seqs.starting_seq = self.starting_seq
        filtered_seqs.name = self.name + '_%s-filtered-%.2f' % ('-'.join(score_types), percentile_cutoff)
        filtered_seqs.seq_dicts = [seq_dict for score, seq_count, seq_dict in filtered_scores]
        filtered_seqs.raw_seqs.extend( [d['sequence'] for d in filtered_seqs.seq_dicts] )
        return filtered_seqs

    def filter_by_seq_position(self, pos, seq, seq_present = True):
        filtered_seqs = Seqs()
        filtered_seqs.alphabet = self.alphabet
        filtered_seqs.alphabet_as_list_shuffled = self.alphabet_as_list_shuffled
        filtered_seqs.raw_seqs = corebio.seq.SeqList( alphabet = self.alphabet )
        filtered_seqs.seq_length = self.seq_length
        filtered_seqs.design_positions = self.design_positions
        filtered_seqs.starting_seq = self.starting_seq
        if seq_present:
            filtered_seqs.name = self.name + '_%d-is%s' % (self.design_positions[pos], seq)
        else:
            filtered_seqs.name = self.name + '_%d-not%s' % (self.design_positions[pos], seq)
        filtered_seqs.seq_dicts = []

        for raw_seq, seq_dict in zip(self.raw_seqs, self.seq_dicts):
            add_seq = False
            if seq_present:
                if raw_seq[pos] == seq:
                    add_seq = True
            else:
                if raw_seq[pos] != seq:
                    add_seq = True

            if add_seq:
                filtered_seqs.raw_seqs.append(raw_seq)
                filtered_seqs.seq_dicts.append(seq_dict)
        return filtered_seqs

    def delta_other(self, other):
        self_matrix = self.get_counts().as_matrix()
        other_matrix = other.get_counts().as_matrix()
        assert( self_matrix.shape == other_matrix.shape )
        return np.subtract(self_matrix, other_matrix)

    def as_matrix(self):
        return self.get_counts().as_matrix()

    def get_counts(self):
        counts = Counts(self.seq_length)
        for raw_seq in self.raw_seqs:
            counts.add_seq(raw_seq)
        return counts

    def __len__(self):
        return len(self.raw_seqs)

    def make_pdf_weblogo(self, pdf_path, prior = None ): # prior could also be = weblogolib.parse_prior('equiprobable',  corebio.seq.unambiguous_protein_alphabet) )
        assert( pdf_path.endswith('.pdf') )
        assert( self.design_positions != None )
        assert( self.starting_seq != None )
        pdf_basedir = os.path.dirname(pdf_path)
        if not os.path.isdir(pdf_basedir):
            os.makedirs(pdf_basedir)
        eps_logo_filename = pdf_path[:-4] + '.eps'
        data = weblogolib.LogoData.from_seqs(self.raw_seqs, prior = prior )
        options = weblogolib.LogoOptions()
        options.show_fineprint = False
        options.xaxis_tic_interval = 1
        options.number_interval = 1
        options.number_fontsize = 3
        options.stacks_per_line = 40
        options.show_errorbars = False
        logo_format = weblogolib.LogoFormat(data, options)
        eps_binary = weblogolib.eps_formatter( data, logo_format )
        eps_str = eps_binary.decode()
        eps_str = replace_logo_numbers(eps_str, self.design_positions, self.starting_seq )
        with open( eps_logo_filename, 'w') as f:
            f.write( eps_str )
        eps_to_pdf(eps_logo_filename)

    def compute_mi(self, use_multiprocessing = True, mi_func='nmi'):
        self.mi = { x : [] for x in xrange(self.seq_length) }

        aln = LoadSeqs( data = { x+1 : y for x,y in enumerate(self.raw_seqs)}, moltype=PROTEIN, aligned=DenseAlignment)
        validate_alignment(aln)


        matrix = coevolve_alignment(coevolve_alignment_functions[mi_func] ,aln)
        for x in xrange( len(matrix) ):
            for y in xrange( len(matrix) ):
                self.mi[x].append( (matrix[x][y], y) )
            self.mi[x].sort( reverse = True )

        # method = coevolve_pair_functions[mi_func]
        # if method == sca_pair:
        #     sca_input_validation(aln)
        # elif method == ancestral_state_pair:
        #     ancestral_states_input_validation(alignment)
        # total_pairs = 0
        # for pos1 in xrange( self.seq_length ):
        #     validate_position(aln, pos1)
        #     for pos2 in xrange( self.seq_length ):
        #         if pos1 <= pos2:
        #             total_pairs += 1

        # if use_multiprocessing:
        #     pool = multiprocessing.Pool( processes = 4 )

        # r = Reporter('calculation mutual information for %s with the "%s" function' % (self.name, mi_func), entries = 'pairs' )
        # r.set_total_count(total_pairs)

        # def callback_helper(val_tup):
        #     pos1, pos2, val = val_tup
        #     for pos, other_pos in [(pos1, pos2), (pos2, pos1)]:
        #         self.mi[pos].append( (val, other_pos) )
        #     r.increment_report()

        # for pos1 in xrange( self.seq_length ):
        #     for pos2 in xrange( self.seq_length ):
        #         if pos1 <= pos2:
        #             if use_multiprocessing:
        #                 pool.apply_async( pool_helper, (pos1, pos2, method, aln), callback = callback_helper )
        #             else:
        #                 callback_helper( pool_helper(pos1, pos2, method, aln) )

        # if use_multiprocessing:
        #     pool.close()
        #     pool.join()
        # r.done()

def pool_helper(pos1, pos2, method, aln):
    return ( pos1, pos2, method(aln, pos1=pos1, pos2=pos2) )

class Counts:
    def __init__(self, num_positions):
        self.alphabet = corebio.seq.unambiguous_protein_alphabet
        self.total_seqs = 0
        self.position_counts = []
        for i in xrange(num_positions):
            self.position_counts.append({})
            for char in self.alphabet:
                self.position_counts[i][char] = 0

    @property
    def num_positions(self):
        return len(self.position_counts)

    @property
    def freqs(self):
        position_freqs = []
        for position_count_dict in self.position_counts:
            d = {}
            for char, count in position_count_dict.iteritems():
                d[char] = float(count) / float(self.total_seqs)
            position_freqs.append(d)
        return position_freqs

    def as_matrix(self):
        matrix = []
        for position_count_dict in self.position_counts:
            matrix_line = []
            for char in sorted([x for x in self.alphabet]):
                matrix_line.append( np.float64(float(position_count_dict[char] / float(self.total_seqs)) ) )
            matrix.append(matrix_line)
        return np.array( matrix )

    def add_seq(self, seq):
        assert( len(seq) == self.num_positions )
        for i, char in enumerate(seq):
            self.position_counts[i][char] += 1
        self.total_seqs += 1

    def display(self):
        print 'Counts object based on %d total sequences:' % self.total_seqs
        for i, freq_d in enumerate(self.freqs):
            freq_str = ' Position %02d - ' % (i+1)
            for char in sorted(freq_d.keys()):
                if freq_d[char] >= 0.01:
                    freq_str += '%s: %.2f, ' % (char, freq_d[char])
            freq_str = freq_str[:-2] # Remove trailing comma and space
            print freq_str

def eps_to_pdf(eps_path):
    assert( eps_path.endswith('.eps') )
    subprocess.check_call( ['epstopdf', eps_path] )
    pdf_path = eps_path[:-4] + '.pdf'
    assert( os.path.isfile(pdf_path) )
    os.remove( eps_path )
    return pdf_path

def parse_rosetta_output( file_path ):
    design_positions_str = 'protocols.coupled_moves: Design Positions:'
    starting_sequence_str = 'Starting Sequence:'
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith(design_positions_str):
                design_positions = [int(x) for x in line[len(design_positions_str):line.find(starting_sequence_str)].strip().split()]
                starting_seq = line[line.find(starting_sequence_str)+len(starting_sequence_str):].strip()
                return design_positions, starting_seq

def replace_logo_numbers(eps_str, sorted_design_positions, ordered_starting_seq ):
    new_eps_str_lines = []
    pos_count = 0
    for line in eps_str.split('\n'):
        if line.endswith('StartStack') and not line.startswith('%'):
            line = line.replace( str(pos_count + 1), '%s%d' % (ordered_starting_seq[pos_count], sorted_design_positions[pos_count]) )
            pos_count += 1
        new_eps_str_lines.append( line)

    return '\n'.join( new_eps_str_lines )

def analyze_single_dir(output_dir, dir_name):
    assert( os.path.isdir(output_dir) )

    stats_files = []
    design_positions = None
    starting_seq = None
    def recursive_search(cwd, stats_files, design_positions, starting_seq):
        for output_file in os.listdir(cwd):
            output_file_path = os.path.join(cwd, output_file)
            if output_file.endswith('.log'):
                if design_positions == None:
                    design_positions, starting_seq = parse_rosetta_output( output_file_path )
                else:
                    this_design_positions, this_starting_seq = parse_rosetta_output( output_file_path )
                    assert( this_design_positions == design_positions )
                    assert( this_starting_seq == starting_seq )
            if output_file.endswith('.stats'):
                stats_files.append(output_file_path)

        for subdir in os.listdir(cwd):
            subdir_path = os.path.join(cwd, subdir)
            if os.path.isdir(subdir_path):
                stats_files, design_positions, starting_seq = recursive_search(subdir_path, stats_files, design_positions, starting_seq)

        return (stats_files, design_positions, starting_seq)

    stats_files, design_positions, starting_seq = recursive_search(output_dir, stats_files, design_positions, starting_seq)

    seqs = Seqs()
    seqs.design_positions = design_positions
    seqs.starting_seq = starting_seq
    seqs.name = dir_name
    seqs.dir_path = os.path.dirname(output_dir)
    for stats_path in stats_files:
        seqs.load_seqs_from_stats(stats_path)
    print '{} sequences loaded from {} stats files'.format(len(seqs), len(stats_files))

    seqs.make_pdf_weblogo( os.path.join(output_dir, '%s-logo.pdf' % seqs.name ), prior = None )

    return seqs

if __name__ == '__main__':
    output_dirs = sys.argv[1:]
    all_seqs = []
    for output_dir in output_dirs:
        assert( os.path.isdir(output_dir) )
        dir_name = os.path.basename(output_dir)
        print 'Analyzing:', dir_name
        all_seqs.append( analyze_single_dir(output_dir, dir_name) )

    topx_to_show = 10

    # Commented out as this tends to give more spurious results
    # wt_seqs = all_seqs[0].make_wt()
    # for seqs in all_seqs:
    #     delta_seqs = Seqs.self_minus_other(seqs, wt_seqs)
    #     delta_seqs.make_pdf_weblogo( os.path.join('plots', '%s_over_wt.pdf' % seqs.name), prior = None )
    #     mutations = seqs.mutations_enriched_over_other(wt_seqs)
    #     print 'Top %d enriched mutations (enriched in %s over wt)' % (topx_to_show, seqs.name)
    #     for enrichment_factor, resi, mut_aa, wt_aa in mutations[:topx_to_show]:
    #         print '%d%s->%s: %.2f' % (resi, wt_aa, mut_aa, enrichment_factor)
    #     print

    # Precompute mutual information
    root_text_output_dir = os.path.join('output', 'text')
    for seqs in all_seqs:
        seqs.compute_mi()
        text_output_dir = os.path.join(root_text_output_dir, seqs.dir_path)
        cov_dump_path = os.path.join(text_output_dir, '%s_covariation.pickle' % (seqs.name) )
        if not os.path.isdir(text_output_dir):
            os.makedirs(text_output_dir)
        with open(cov_dump_path, 'w') as f:
            pickle.dump(seqs.mi, f)


    for i, seqs_pair0 in enumerate(all_seqs):
        for j, seqs_pair1 in enumerate(all_seqs):
            if i != j:
                # scores_to_filter = ['total_score'] # ['angle_constraint', 'atom_pair_constraint']
                # percentile_filter = 0.7
                # filtered0 = seqs_pair0.filter_by_score_percentile(scores_to_filter, percentile_filter)
                # filtered1 = seqs_pair1.filter_by_score_percentile(scores_to_filter, percentile_filter)
                seq_pairs = [
                    (seqs_pair0, seqs_pair1),
                    # (filtered0, filtered1)
                ]
                for seqs0, seqs1 in seq_pairs:
                    delta_seqs = Seqs.self_minus_other( seqs0, seqs1 )
                    delta_seqs.make_pdf_weblogo( os.path.join('plots', '%s.pdf' % delta_seqs.name), prior = None )
                    mutations = seqs0.mutations_enriched_over_other(seqs1)
                    print 'Top %d enriched mutations (enriched in %s over %s)' % (topx_to_show, seqs0.name, seqs1.name)
                    printed_lines = 0

                    assert( seqs0.dir_path == seqs1.dir_path )
                    with open(os.path.join(text_output_dir, '%s_over_%s.tsv' % (seqs0.name, seqs1.name) ), 'w') as f:
                        f.write('enrichment_factor\twt_aa\tresnum\tmut_aa\t%s_freq\t%s_freq\n' % (seqs0.name, seqs1.name) )

                        # Get top mutations
                        top_enriched_mutations_index = []
                        top_enriched_mutations_restype = []
                        for enrichment_factor, resi, mut_aa, wt_aa, freq0, freq1, raw_index in mutations:
                            if mut_aa != wt_aa:
                                top_enriched_mutations_index.append(raw_index)
                                top_enriched_mutations_restype.append(mut_aa)


                        for enrichment_factor, resi, mut_aa, wt_aa, freq0, freq1, raw_index in mutations:
                            if mut_aa == wt_aa:
                                continue
                            mutation_description = '%d%s->%s: %.2f (%s freq: %.2f; %s freq: %.2f)' % (resi, wt_aa, mut_aa, enrichment_factor, seqs0.name, freq0, seqs1.name, freq1)
                            covariation_cache = {}
                            for cov_val, y in seqs0.mi[raw_index]:
                                if y !=raw_index:
                                    covariation_cache[seqs0.design_positions[y]] = cov_val

                            # Look for enrichment at other postions in sequences with this mutation
                            mut_present = seqs0.filter_by_seq_position(raw_index, mut_aa, seq_present = True)
                            mut_absent = seqs0.filter_by_seq_position(raw_index, mut_aa, seq_present = False)
                            mut_present_lines = []
                            min_present_count = min( len(mut_present), len(mut_absent) )
                            if min_present_count > 0:
                                enriched_for_mut = mut_present.mutations_enriched_over_other(mut_absent)
                                for mut_present_enrichment_factor, mut_present_resi, mut_present_mut_aa, mut_present_wt_aa, mut_present_freq0, mut_present_freq1, mut_present_raw_index in enriched_for_mut:
                                    if mut_present_mut_aa != mut_present_wt_aa and not (mut_present_mut_aa == mut_aa and mut_present_wt_aa == wt_aa) and len(mut_present_lines) < 10 and abs(mut_present_enrichment_factor) >= 0.05:
                                        if mut_present_resi in covariation_cache:
                                            mut_present_mutation_description = '%03d%s->%s: %.2f (%s freq: %.2f, n=%d; %s freq: %.2f, n=%d) covariation=%.2f' % (mut_present_resi, mut_present_wt_aa, mut_present_mut_aa, mut_present_enrichment_factor, mut_present.name, mut_present_freq0, len(mut_present), mut_absent.name, mut_present_freq1, len(mut_absent), covariation_cache[mut_present_resi])
                                        else:
                                            mut_present_mutation_description = '%03d%s->%s: %.2f (%s freq: %.2f, n=%d; %s freq: %.2f, n=%d) covariation=N/A' % (mut_present_resi, mut_present_wt_aa, mut_present_mut_aa, mut_present_enrichment_factor, mut_present.name, mut_present_freq0, len(mut_present), mut_absent.name, mut_present_freq1, len(mut_absent) )
                                        mut_present_lines.append(mut_present_mutation_description)
                                if len(mut_present_lines) > 1:
                                    mut_present_lines.insert(0, 'Mutations enriched in %s over %s when %d%s->%s is present' % (seqs0.name, seqs1.name, resi, wt_aa, mut_aa) )

                            if abs(enrichment_factor) >= 0.01:
                                f.write( '%.2f\t%s\t%d\t%s\t%.2f\t%.2f\n' % (enrichment_factor, wt_aa, resi, mut_aa, freq0, freq1) )
                            if printed_lines < topx_to_show:
                                print mutation_description
                                for line in mut_present_lines:
                                    print line
                                print
                                printed_lines += 1
                    print
