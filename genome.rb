#
# genome.rb
# A tool for inspecting FASTA-format genomes.
#
#
#
# alexander.s.farley@gmail.com
# May 24, 2020
#

require './gene.rb'
require 'colorize'

INPUT_FILE_PATH = ARGV[0]
input_file_lines = IO.readlines(INPUT_FILE_PATH)
input_file_header = input_file_lines[0]
input_file_genome_lines = input_file_lines[1..-1]
genome = input_file_genome_lines.join("").gsub(/[^0-9a-z ]/i, '')
genome_scan_variations = genome_to_scan_variations(genome)
genes = Set.new # Declaring genes as a Set (rather than an array) means that we don't have to search for inclusion when appending genes from each scan.

genome_scan_variations.each do |genome_scan|
	# Split scan into triplets for comparison with start and end-codons. Discard elements containing less than three bases.
	triplets = genome_scan.to_triplets
	
	# Scan triplets for matches to [START_CODON ... END_CODON], including nested ORFs with multiple start-codons (but not multiple end-codons).
	genes_in_scan = triplets.loose_start_hard_end(START_CODON,END_CODONS)
	genes_in_scan.each{ |gene| genes << gene.to_amino if gene.length >= MIN_ORF_LENGTH }
end

genes.each do  |gene|
	puts "#{gene}\r\n"
	puts "Length:#{gene.length}\r\n\r\n".yellow
end