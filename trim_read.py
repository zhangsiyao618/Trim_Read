import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import extend_adapter


class FASTAQ:
    def __init__(self, input_file):
        self.filename = input_file
        self.sequence_lines = list()
        self.scores = list()
        self.records = list()
        self.score_encoding = ''
        self.phred_scores = list()
        self.info = list()
        self.length = list()
    
    # Extract information from the FASTAQ file
    def extract_information(self):
        with open(self.filename, 'r') as file:
            lines = file.readlines()
        file.close() 

        for line in lines[1::4]:
            line = line.strip('\n')
            self.sequence_lines.append(line)
            self.length.append(len(line))
        for line in lines[0::4]:
            line = line.strip('\n')
            self.info.append(line)
        for record in lines[3::4]:
            score = list()
            record = record.strip('\n')
            self.records.append(record)
            for char in record:
                score.append(ord(char))
            self.scores.append(score)
        file.close()
    
    # Determine the score encoding system based on min/max values and convert the score encodings to Phred33 or Phred64 or Solexa64 according to the score encoding
    def decode(self):
        # Combine the scores into a single list 
        combined_scores = list()
        for score in self.scores:
            combined_scores.extend(score)
        min_score = min(combined_scores)
        max_score = max(combined_scores)

        # Determine the score encoding system based on min/max values
        if max_score <= 126 and min_score < 59:
            self.score_encoding = 'Phred33'
        elif max_score > 73 and min_score >= 64:
            self.score_encoding = 'Phred64'
        elif min_score >= 59 and min_score < 64 and max_score > 73:
            self.score_encoding = 'Solexa64'
        else:
            self.score_encoding = ('Unknown score encoding', min_score, max_score)

        # Convert the score encodings to Phred33 or Phred64 or Solexa64 according to the score encoding
        self.phred_scores = list()
        if self.score_encoding == 'Phred33':
            for score in self.scores:
                phred_score = list()
                for i in score:
                    phred_score.append(i-33)
                self.phred_scores.append(phred_score)
        elif self.score_encoding == 'Phred64' or self.score_encoding == 'Solexa64':
            for score in self.scores:
                phred_score = list()
                for i in score:
                    phred_score.append(i-64)
                self.phred_scores.append(phred_score)
        else:
            print('Unknown score encoding')
        print('Encoding system :', self.score_encoding)


class Quality_Assessment:
    def __init__(self, sequences, records, phred_scores, length):
        self.sequence_lines = sequences
        self.records = records
        self.phred_scores = phred_scores
        self.length = length

    # Plot the distribution of sequence lengths
    def length_distribution(self):
        plt.hist(self.length, bins=20, color='steelblue', edgecolor='black')
        plt.xlabel('Length')
        plt.ylabel('Frequency')
        plt.title('Sequence Length Distribution')
        plt.show()
    
    # Produce an assessment table of quality statistics 
    def assess_quality(self):
        # Create a dataframe to store the Phred scores
        self.phred_df = pd.DataFrame({'Read ID': range(1, len(self.phred_scores)+1)})
        # Add the Phred scores to the dataframe
        self.phred_df['Phred Scores'] = pd.Series(self.phred_scores)
        self.phred_df['Mean'] = self.phred_df['Phred Scores'].apply(lambda x: np.mean(x) if len(x) > 0 else np.nan)
        self.phred_df['Std'] = self.phred_df['Phred Scores'].apply(lambda x: np.std(x) if len(x) > 0 else np.nan)
        self.phred_df['Median'] = self.phred_df['Phred Scores'].apply(lambda x: np.median(x) if len(x) > 0 else np.nan)
        self.phred_df['Min'] = self.phred_df['Phred Scores'].apply(lambda x: np.min(x) if len(x) > 0 else np.nan)
        self.phred_df['Max'] = self.phred_df['Phred Scores'].apply(lambda x: np.max(x) if len(x) > 0 else np.nan)

        phred_qualities = list()
        single_quality_counts = list()
        self.whole_quality_count = dict()
        for phred_score in self.phred_scores:
            # Grade the Phred scores
            phred_quality = pd.cut(phred_score, bins=[-1, 10, 20, 30, 100], labels=['Poor', 'Medium', 'Good', 'Excellent'])
      
            # Count the frequency of quality grades for each sequence
            single_quality_count = list()
            single_quality_count.append(list(phred_quality).count('Poor'))
            single_quality_count.append(list(phred_quality).count('Medium'))
            single_quality_count.append(list(phred_quality).count('Good'))
            single_quality_count.append(list(phred_quality).count('Excellent'))

            # Count the frequency of quality grades for all reads
            self.whole_quality_count['Low'] = self.whole_quality_count.get('Low', 0) + list(phred_quality).count('Poor')
            self.whole_quality_count['Medium'] = self.whole_quality_count.get('Medium', 0) + list(phred_quality).count('Medium')
            self.whole_quality_count['Good'] = self.whole_quality_count.get('Good', 0) + list(phred_quality).count('Good')
            self.whole_quality_count['Excellent'] = self.whole_quality_count.get('Excellent', 0) + list(phred_quality).count('Excellent')
        
            phred_qualities.append(phred_quality) 
            single_quality_counts.append(single_quality_count)

        self.phred_df['Quality'] = pd.Series(phred_qualities)
        self.phred_df['Quality_cnt'] = pd.Series(single_quality_counts)

    # Calculate gc_content(per site) of the original sequences and draw charts of GC Content Distribution
    def GC_content(self):
        gc_contents = []
        # Transpose the sequence matrix and calculate the GC content per site of the reads
        sequence_lines = [x for x in self.sequence_lines if x != '']
        transposed_sequences = ["".join(sequence_line) for sequence_line in zip(*sequence_lines)]
        
        for sequence_line in transposed_sequences:
            sequence = sequence_line.strip()
            gc_count = sequence.count("G") + sequence.count("C")
            gc_content = gc_count / len(sequence) * 100
            gc_contents.append(gc_content)
        self.gc_count = np.histogram(gc_contents, bins=20, range=(0, 100))[0]
        
        x = np.linspace(0, 100, 20)
        y = self.gc_count
        fig, ax = plt.subplots()
        # Draw a histogram and a line chart of GC Content Distribution
        ax.bar(x, y, width=5, edgecolor='black', alpha=0.7, label='Histogram')
        ax.plot(x, y, marker='o', linestyle='-', color='red', label='Line')
        ax.set_xlabel('Per site GC Content (%)')
        ax.set_ylabel('Count')
        ax.set_title('GC Content Distribution')
        ax.legend()
        plt.show()

    # Calculate the phred scrores of all sites and draw a scatter plot
    def site_scatter(self):
        # Generate the data of the abscissa (site) and the ordinate (Phred score)
        combined_scores = list()
        for score in self.phred_scores:
            if len(score) == 0:
                continue
            combined_scores.extend(score)
        positions = range(1, len(combined_scores) + 1)
     
        # Draw a scatter plot
        plt.scatter(positions, combined_scores, s=10, c='green', alpha=0.5)
        plt.title("Phred Score Distribution")
        plt.xlabel("Position")
        plt.ylabel("Phred Score")
        plt.show()

    def seq_quality_per_site(self):
        phred_scores = [x for x in self.phred_scores if x != []]
        # Transpose the phred score matrix and calculate the phred score per site of the reads
        transposed_scores = [list(row) for row in zip(*phred_scores)]
        # Calculate the mean phred score of bases each site
        mean_phred_scores = [sum(column) / len(column) for column in transposed_scores]

        # Draw a scatter plot
        plt.scatter(range(len(mean_phred_scores)), mean_phred_scores, s=10, c='blue')
        plt.xlabel('Position')
        plt.ylabel('Mean Quality Score')
        plt.title('Sequence Quality Per Site ')
        plt.show()


class TrimRead:
    def __init__(self, reads, phred_scores, records, adapter_list, window_size, threshold1, threshold2):
        self.reads = reads
        self.phred_scores = phred_scores
        self.records = records
        self.window_size = window_size
        self.threshold1 = threshold1
        self.threshold2 = threshold2
        self.adapter_list = adapter_list
    
    def trim(self):
        self.cut_reads, self.cut_phred_scores, self.cut_records = self.cut_adapter(self.reads, self.phred_scores, self.records, self.adapter_list, self.threshold1)
        if len(self.reads) < 500:
            trimmed_5_reads, trimmed_5_phred_scores, trimmed_5_records = self.base_by_base5(self.cut_reads, self.cut_phred_scores, self.cut_records, self.threshold2)
            self.trimmed_reads, self.trimmed_phred_scores, self.trimmed_records = self.base_by_base3(trimmed_5_reads, trimmed_5_phred_scores, trimmed_5_records, self.threshold2)
        else:
            trimmed_5_reads, trimmed_5_phred_scores, trimmed_5_records = self.window_by_window5(self.cut_reads, self.cut_phred_scores, self.cut_records, self.window_size, self.threshold2)
            self.trimmed_reads, self.trimmed_phred_scores, self.trimmed_records = self.window_by_window3(trimmed_5_reads, trimmed_5_phred_scores, trimmed_5_records, self.window_size, self.threshold2)

    @staticmethod
    # Cut the adapter from the 5' end of read
    def cut_adapter(reads, phred_scores, records, adapter_list, threshold):
        cut_reads = []
        cut_phred_scores = []
        cut_records = []
    
        for read, phred_score, record in zip(reads, phred_scores, records):
            error_rate_list = []
            alignment_list = []
            match_len_list = []
            for adapter in adapter_list:
                error_rate, alignment, match_len = TrimRead.find_match(read[0:(len(adapter)+5)], adapter)
                if match_len >= len(adapter) * (1 - threshold) and error_rate < threshold:
                    error_rate_list.append(error_rate)
                    alignment_list.append(alignment)
                    match_len_list.append(match_len)
                else:
                    continue    

            if len(error_rate_list) != 0:  # if there is a match satisfying the threshold set on length and error rate, cut the match from read
                lowest_error_rate = min(error_rate_list)
                best_match = alignment_list[error_rate_list.index(lowest_error_rate)]
                cut_length = match_len_list[error_rate_list.index(lowest_error_rate)]
                start = read.index((best_match[0]).replace('-', ''))
                if start > 0:
                    start = 0
                cut_read = read[:start] + read[(start + cut_length):]
                cut_phred_score = phred_score[:start] + phred_score[(start + cut_length):]
                cut_record = record[:start] + record[(start + cut_length):]
                cut_reads.append(''.join(cut_read))
                cut_phred_scores.append(cut_phred_score)
                cut_records.append(''.join(cut_record))
            else:
                cut_reads.append(read)
                cut_phred_scores.append(phred_score)
                cut_records.append(record)
        
        return cut_reads, cut_phred_scores, cut_records

    @staticmethod
    # Find the best match between adapter and read and return the result
    def find_match(seq1, seq2, match_score=1, gap_penalty=-1, mismatch_penalty=-1):
        # Use the Smith-Waterman algorithm to find the best match between adapter and read
        m = len(seq1)
        n = len(seq2)
        i_max, j_max = 0, 0
        max_score = 0
        matrix = [[0] * (n + 1) for _ in range(m + 1)]
        for i in range(1, m+1):
            for j in range(1, n+1):
                if seq1[i-1] == seq2[j-1]:
                    match = matrix[i-1][j-1] + match_score
                else:
                    match = matrix[i-1][j-1] + mismatch_penalty
                delete = matrix[i-1][j] + gap_penalty 
                insert = matrix[i][j-1] + gap_penalty
                matrix[i][j] = max(match, delete, insert, 0)
                if matrix[i][j] > max_score:
                    max_score = matrix[i][j]
        for i in range(1, m+1):
            for j in range(1, n+1):
                if matrix[i][j] == max_score:
                    i_max = i
                    j_max = j

        # Traceback to find the best match
        match_length, mismatch_length = 0, 0
        i, j = i_max, j_max
        align1, align2 = '', ''
        while matrix[i][j] > 0:
            score = matrix[i][j]
            up = matrix[i-1][j]
            left = matrix[i][j-1]
            diag = matrix[i-1][j-1]
            if score == diag + match_score or score == diag + mismatch_penalty:
                align1 = seq1[i-1] + align1
                align2 = seq2[j-1] + align2
                i -= 1
                j -= 1
                if score == diag + mismatch_penalty :
                    mismatch_length += 1
                else:
                    match_length += 1
            elif score == up + gap_penalty:
                align1 = seq1[i-1] + align1
                align2 = "-" + align2
                i -= 1
            elif score == left + gap_penalty:
                align1 = "-" + align1
                align2 = seq2[j-1] + align2
                j -= 1
        alignment = ["".join(align1), "".join(align2)]
        match_length = len(align1)-str.count(align1, '-')
        error_rate = mismatch_length / match_length if match_length != 0 else 1
        return error_rate, alignment, match_length
    
    @staticmethod
    # Trim the 5' end of the sequence in a base-by-base way
    def base_by_base5(reads, phred_scores, records, threshold):
        trimmed_read = []
        trimmed_phred_score = []
        trimmed_records = []
        for read, phred_score, record in zip(reads, phred_scores, records):
            # If the sequence is already empty, return empty string
            if len(read) == 0:
                trimmed_read.append('')
                trimmed_phred_score.append([])
                trimmed_records.append('')
                continue

            # Calculate differences between quality scores and threshold
            diffs = [q - threshold for q in phred_score]

            # Calculate cumulative sum of differences
            cumsum = [sum(diffs[:i+1]) for i in range(len(diffs))]

            # Find minimum cumulative sum index
            idx = cumsum.index(min(cumsum))

            # Trim sequence and quality lists
            trimmed_read.append(read[idx+1:])
            trimmed_phred_score.append(phred_score[idx+1:])
            trimmed_records.append(record[idx+1:])
        return trimmed_read, trimmed_phred_score, trimmed_records

    @staticmethod
    # Trim the 3' end of the sequence in a base-by-base way
    def base_by_base3(reads, phred_scores,records, threshold):
        # Convert quality string to list of integers
        trimmed_read = []
        trimmed_phred_score = []
        trimmed_records = []
        for read, phred_score, record in zip(reads, phred_scores, records):
            # If the sequence is already empty, return empty string
            if len(read) == 0:
                trimmed_read.append('')
                trimmed_phred_score.append([])
                trimmed_records.append('')
                continue
        
            # Calculate differences between quality scores and threshold 计算与阈值差值
            diffs = [q - threshold for q in reversed(phred_score)]  
            # Calculate cumulative sum of differences 计算差值累计
            cumsum = list(reversed([sum(diffs[i:]) for i in range(len(diffs))]))

            # Find minimum cumulative sum index 找到累计差值最小处
            idx = cumsum.index(min(cumsum))

            # Trim sequence and quality lists 在累计差值为最小处截断
            trimmed_read.append(read[:len(read)-idx - 1])
            trimmed_phred_score.append(phred_score[:len(phred_score)-idx - 1])
            trimmed_records.append(record[:len(record)-idx - 1])
        return trimmed_read, trimmed_phred_score, trimmed_records
    
    @staticmethod
    # Trim the 5' end of the sequence in a base-by-base way
    def window_by_window5(reads, phred_scores, records, window_size, threshold):
        trimmed_read = []
        trimmed_phred_score = []
        trimmed_records = []
        for read, phred_score, record in zip(reads, phred_scores, records):
            # Set initial values for variables
            start = 0
            end = len(read)
            mean_qual = sum(phred_score[start:end]) / float(end - start)

            # Move sliding window to the left side of the sequence
            while end - start > window_size and mean_qual < threshold:
                # Calculate quality for the previous window
                next_end = start + window_size
                next_start = start
                next_mean_qual = sum(phred_score[next_start:next_end]) / float(window_size)

                # If previous window quality is above threshold, stop trimming
                if next_mean_qual >= threshold:
                    break

                # Trim current window
                start = next_end
                mean_qual = next_mean_qual
        
            # If the left read is too short, then we will not keep it
            if end - start < len(read) * 0.5 or end - start < 20:
                trimmed_read.append('')
                trimmed_phred_score.append([])
                trimmed_records.append('')
            else:
                # Return trimmed sequence and quality lists
                trimmed_read.append(read[start:])
                trimmed_phred_score.append(phred_score[start:])
                trimmed_records.append(record[start:])
        return trimmed_read, trimmed_phred_score, trimmed_records
    
    @staticmethod
    # Trim the 3' end of the sequence in a window-by-window way
    def window_by_window3(reads, phred_scores, records, window_size, threshold):
        trimmed_read = []
        trimmed_phred_score = []
        trimmed_records = []
        for read, phred_score, record in zip(reads, phred_scores, records):
            # Set initial values for variables
            start = 0
            end = len(read)
            mean_qual = sum(phred_score[start:end]) / float(end - start) if end != 0 else 0

            # Move sliding window to the right side of the sequence
            while end - start > window_size and mean_qual < threshold:
                # Calculate quality for the next window
                next_start = end - window_size
                next_end = end
                next_mean_qual = sum(phred_score[next_start:next_end]) / float(window_size)

                # If next window quality is above threshold, stop trimming
                if next_mean_qual >= threshold:
                    break

                # Trim current window
                end = next_start
                mean_qual = next_mean_qual
        
            # If the left read is too short, then we will not keep it
            if end - start < len(read) * 0.5 or end - start < 30:
                trimmed_read.append('')
                trimmed_phred_score.append([])
                trimmed_records.append('')
            else:
                # Return trimmed sequence and quality lists
                trimmed_read.append(read[:end])
                trimmed_phred_score.append(phred_score[:end])
                trimmed_records.append(record[:end])

        return trimmed_read, trimmed_phred_score, trimmed_records


# Write the cut reads to a new fastq file
def write_cut_reads_fastq(cut_reads, cut_records, info, output_file):
    with open(output_file, 'w') as f:
        i = 0
        for read in cut_reads:
            if read == '':
                i += 1
                continue
            f.write(info[i])
            f.write('\n')
            f.write(read)
            f.write('\n')
            f.write('+')
            f.write('\n')
            f.write(cut_records[i])
            f.write('\n')
            i += 1

# Write the cut reads to a new fasta file
def write_cut_reads_fasta(cut_reads, output_file):
    with open(output_file,'w') as f:
        i = 1
        for read in cut_reads:
            f.write('>Seq'+str(i)+'\n')
            if read == '':
                f.write('Abandoned for low quality.'+'\n')
            else:
                f.write(read+'\n')
            i += 1

# Read the adapter sequences from a file
def read_file_lines(filename):
    lines = []
    with open(filename, 'r') as file:
        for line in file:
            lines.append(line.strip())
    return lines

# Create a parser object
parser = argparse.ArgumentParser()
parser.add_argument('input_filepath', type=str, help='The untreated input file.')  
parser.add_argument('--fastq', '-q', action='store_true', help='The output format is chosen to be fastq.')
parser.add_argument('--fasta', '-a', action='store_true', help='The output format is chosen to be fasta.')
parser.add_argument('--window_size', '-w', type=int, default=5, help='The sliding window size.')
parser.add_argument('--threshold1', '-t1', type=float, default=0.2, help='Threshold for error rate.')
parser.add_argument('--threshold2', '-t2', type=float, default=10, help='Threshold for quality scores.')
parser.add_argument('--adapter5_list', '-l', type=str, default='adapter5_list.txt', help='file:adapter5_list')

# Parse the arguments from standard input
args = parser.parse_args()

extend_adapter.reverse_complement_each_line(args.adapter5_list)
adapter_list = read_file_lines(args.adapter5_list)

input_file = args.input_filepath
object = FASTAQ(input_file)
object.extract_information()
object.decode()

# Quality accessment before trimming
result_original = Quality_Assessment(object.sequence_lines, object.records, object.phred_scores, object.length)
result_original.length_distribution()
result_original.GC_content()    
result_original.assess_quality()
print('phred_dataframe\n', result_original.phred_df)
result_original.site_scatter()
result_original.seq_quality_per_site()

# Trim the reads
trim_read = TrimRead(object.sequence_lines, object.phred_scores, object.records, adapter_list, args.window_size, args.threshold1, args.threshold2)
trim_read.trim()

# Quality accessment after trimming
trimmed_length = []
for read in trim_read.trimmed_reads:
    trimmed_length.append(len(read))
result_after = Quality_Assessment(trim_read.trimmed_reads, trim_read.trimmed_records, trim_read.trimmed_phred_scores, trimmed_length)
result_after.length_distribution()
result_after.GC_content()    
result_after.assess_quality()
print('phred_dataframe\n', result_after.phred_df)
result_after.site_scatter()
result_after.seq_quality_per_site()

# Write the cut reads to a new fastq/fasta file
output_file = 'cut_reads.fastq'
if args.fasta:
    write_cut_reads_fastq(trim_read.trimmed_reads, trim_read.trimmed_records, object.info, output_file)
if args.fastq:
    write_cut_reads_fasta(trim_read.trimmed_reads, output_file)