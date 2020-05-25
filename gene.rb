#
# gene.rb
# A library of accessory functions for inspecting FASTA-format genomes.
#
#
#
# alexander.s.farley@gmail.com
# May 24, 2020
#

require 'set' 
#require 'colorize'

START_CODONS = ["ATG", "GTC"]
END_CODONS = ["TAA","TGA","TAG"]
MIN_ORF_LENGTH = 96

GENETIC_CODE = { "UUU" => 'F', \
"UUC" => 'F', \
"UUA" => 'L', \
"UUG" => 'L', \
"CUU" => 'L', \
"CUC" => 'L', \
"CUA" => 'L', \
"CUG" => 'L', \
"AUU" => 'I', \
"AUC" => 'I', \
"AUA" => 'I', \
"AUG" => 'M', \
"GUU" => 'V', \
"GUC" => 'V', \
"GUA" => 'V', \
"GUG" => 'V', \
"UCU" => 'S', \
"UCC" => 'S', \
"UCA" => 'S', \
"UCG" => 'S', \
"CCU" => 'P', \
"CCC" => 'P', \
"CCA" => 'P', \
"CCG" => 'P', \
"ACU" => 'T', \
"ACC" => 'T', \
"ACA" => 'T', \
"ACG" => 'T', \
"GCU" => 'A', \
"GCC" => 'A', \
"GCA" => 'A', \
"GCG" => 'A', \
"UAU" => 'Y', \
"UAC" => 'Y', \
"UAA" => '', \
"UAG" => '', \
"CAU" => 'H', \
"CAC" => 'H', \
"CAA" => 'Q', \
"CAG" => 'Q', \
"AAU" => 'N', \
"AAC" => 'N', \
"AAA" => 'K', \
"AAG" => 'K', \
"GAU" => 'D', \
"GAC" => 'D', \
"GAA" => 'E', \
"GAG" => 'E', \
"UGU" => 'C', \
"UGC" => 'C', \
"UGA" => '',  \
"UGG" => 'W', \
"CGU" => 'R', \
"CGC" => 'R', \
"CGA" => 'R', \
"CGG" => 'R', \
"AGU" => 'S', \
"AGC" => 'S', \
"AGA" => 'R', \
"AGG" => 'R', \
"GGU" => 'G', \
"GGC" => 'G', \
"GGA" => 'G', \
"GGG" => 'G' }

# A TripletArray is an array of base triplets:
# [ "ATG", "AAA", "TAA"]
class TripletArray < Array
	# Get all matching subsets of TripletArrays contained in this TripletArray.
	def loose_start_hard_end(start_set,end_set)
		self.inject([]) { |accum,elem| 
			accum.each{ |partial| partial << elem if not end_set.include? partial.last }
			accum << TripletArray.new([elem]) if start_set.include? elem
			accum 
		}
	end
	
	def to_protein
		amino_bases = self.map{ |triplet| GENETIC_CODE[triplet.to_RNA] }.join
		(amino_bases[0] == 'V') ? amino_bases[1..] : amino_bases
	end
end

class String
	def to_triplets
		t = TripletArray.new
		self.chars.each_slice(3).select{ |elem| elem.length == 3}.map{ |a| a.join }.each{ |elem| t << elem }
		return t
	end
	
	def to_RNA
		self.gsub("T","U")
	end

	# Generate six different scans of a genome, using 0/1/2 offsets and forward/reverse directions.
	def to_scan_variations
		forward_duplicated_0 = self + self
		forward_duplicated_1 = forward_duplicated_0[1..]
		forward_duplicated_2 = forward_duplicated_0[2..]
		reverse_duplicated_0 = forward_duplicated_0.reverse.to_DNA_complement
		reverse_duplicated_1 = forward_duplicated_1.reverse.to_DNA_complement
		reverse_duplicated_2 = forward_duplicated_2.reverse.to_DNA_complement

		scan_variations = [forward_duplicated_0, \
		forward_duplicated_1, forward_duplicated_2, \
		reverse_duplicated_0, reverse_duplicated_1, \
		reverse_duplicated_2]
	end

	def to_DNA_complement
		complements = {
			"A" => "T",
			"T" => "A",
			"G" => "C",
			"C" => "G"
		}
		self.gsub(/\w/,complements)
	end
end

# Example usage: load a file from the first command-line argument and parse ORFs.
if __FILE__ == $0
	INPUT_FILE_PATH = ARGV[0]
	input_file_lines = IO.readlines(INPUT_FILE_PATH)
	#input_file_header = input_file_lines[0]
	input_file_genome_lines = input_file_lines[1..-1]
	genome = input_file_genome_lines.join("").gsub(/[^0-9a-z ]/i, '')
	scan_variations = genome.to_scan_variations
	genes = Set.new # Declaring genes as a Set (rather than an array) means that we don't have to search for inclusion when appending genes from each scan.

	scan_variations.each do |scan|
		# Split scan into triplets for comparison with start and end-codons. Discard elements containing less than three bases.
		triplets = scan.to_triplets
		#puts "All possible amino-acid chains:"
		#puts "#{triplets.to_protein}\r\n\r\n"

		#puts "Triplets: \r\n#{triplets}\r\n"
		# Scan triplets for matches to [START_CODON ... END_CODON], including nested ORFs with multiple start-codons (but not multiple end-codons).
		genes_in_scan = triplets.loose_start_hard_end(START_CODONS,END_CODONS)
		# Take ORFs above minimum length and covert to proteins (amino-acid chains).
		genes_in_scan.each{ |gene| genes << gene.to_protein if gene.length >= MIN_ORF_LENGTH }
	end

	genes.each do  |gene|
		puts "#{gene}"
		puts "Length:#{gene.length}\r\n"
	end
end