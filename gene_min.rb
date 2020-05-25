require 'set' 
START_CODONS = ["ATG", "GTC"]
END_CODONS = ["TAA","TGA","TAG"]
MIN_ORF_LENGTH = 96
GENETIC_CODE = { "UUU" => 'F', "UUC" => 'F', "UUA" => 'L', "UUG" => 'L', "CUU" => 'L', "CUC" => 'L', "CUA" => 'L', "CUG" => 'L', "AUU" => 'I', "AUC" => 'I', "AUA" => 'I', "AUG" => 'M', "GUU" => 'V', "GUC" => 'V', "GUA" => 'V', "GUG" => 'V', "UCU" => 'S', "UCC" => 'S', "UCA" => 'S', "UCG" => 'S', "CCU" => 'P', "CCC" => 'P', "CCA" => 'P', "CCG" => 'P', "ACU" => 'T', "ACC" => 'T', "ACA" => 'T', "ACG" => 'T', "GCU" => 'A', "GCC" => 'A', "GCA" => 'A', "GCG" => 'A', "UAU" => 'Y', "UAC" => 'Y', "UAA" => '', "UAG" => '', "CAU" => 'H', "CAC" => 'H', "CAA" => 'Q', "CAG" => 'Q', "AAU" => 'N', "AAC" => 'N', "AAA" => 'K', "AAG" => 'K', "GAU" => 'D', "GAC" => 'D', "GAA" => 'E', "GAG" => 'E', "UGU" => 'C', "UGC" => 'C', "UGA" => '',  "UGG" => 'W', "CGU" => 'R', "CGC" => 'R', "CGA" => 'R', "CGG" => 'R', "AGU" => 'S', "AGC" => 'S', "AGA" => 'R', "AGG" => 'R', "GGU" => 'G', "GGC" => 'G', "GGA" => 'G', "GGG" => 'G' }
class TripletArray < Array
	def loose_start_hard_end(start_set,end_set)
		self.inject([]) { |accum,elem| 
			accum.each{ |partial| partial << elem if not end_set.include? partial.last }
			accum << TripletArray.new([elem]) if start_set.include? elem
			accum 
		}
	end
	def to_protein
		amino_bases = self.map{ |triplet| GENETIC_CODE[triplet.gsub("T","U")] }.join
		(amino_bases[0] == 'V') ? amino_bases[1..] : amino_bases
	end
end
class String
	def to_triplets
		t = TripletArray.new
		self.chars.each_slice(3).select{ |elem| elem.length == 3}.map{ |a| a.join }.each{ |elem| t << elem }
		return t
	end
	def to_scan_variations
		forward_0 = self
		forward_1 = self[1..]
		forward_2 = self[2..]
		reverse_0 = self.reverse.to_DNA_complement
		reverse_1 = self.reverse.to_DNA_complement[1..]
		reverse_2 = self.reverse.to_DNA_complement[2..]
		scan_variations = [forward_0,forward_1, forward_2,reverse_0, reverse_1,reverse_2]
	end
	def to_DNA_complement
		complements = {"A" => "T", "T" => "A", "G" => "C", "C" => "G" }
		self.gsub(/\w/,complements)
	end
end
if __FILE__ == $0
	genome = IO.readlines(ARGV[0])[1..-1].join("").gsub(/[^0-9a-z ]/i, '')
	scan_variations = genome.to_scan_variations
	genes = Set.new
	scan_variations.each do |scan|
		genes_in_scan = scan.to_triplets.loose_start_hard_end(START_CODONS,END_CODONS)
		genes_in_scan.each{ |gene| genes << gene.to_protein if gene.length >= MIN_ORF_LENGTH }
	end
	genes.each { |gene| puts "#{gene}\r\n\r\n" }
end