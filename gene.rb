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

START_CODON = "ATG"
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
	def loose_start_hard_end(start_element,end_set)
		self.inject([]) { |accum,elem| 
			accum.each{ |partial| 
				partial << elem if not end_set.include? partial.last
			}
			accum << TripletArray.new([elem]) if elem == start_element
			accum
		}
	end
	
	def to_amino
		self.map{ |triplet| GENETIC_CODE[triplet.to_RNA] }.join
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
end

# Treat the genome as 'circular'; allow matches from one end to the other. Also allow matches in reverse.
def genome_to_scan_variations(g)
	forward_duplicated_0 = g + g
	forward_duplicated_1 = forward_duplicated_0[1..]
	forward_duplicated_2 = forward_duplicated_0[2..]
	reverse_duplicated_0 = forward_duplicated_0.reverse
	reverse_duplicated_1 = forward_duplicated_1.reverse
	reverse_duplicated_2 = forward_duplicated_2.reverse

	scan_variations = [forward_duplicated_0, \
	forward_duplicated_1, forward_duplicated_2, \
	reverse_duplicated_0, reverse_duplicated_1, \
	reverse_duplicated_2]
end