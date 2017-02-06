mrnas = {}
cds = {}
genes = {}

DEFAULT = 0
MRNA = 1
CDS = 2
GENE = 3

state = DEFAULT

path = "/home/eve/Documents/biology/load_data/Anolis_carolinensis/Anolis_1.gbk"
file = File.open(path, "r")

while (line = file.gets)
	prefix = line[0..20].strip
	#break if prefix == "ORIGIN"

	if prefix.end_with?("RNA")
		state = MRNA
	elsif prefix == "CDS"
		state = CDS
	elsif prefix == "gene"
		state = GENE
	elsif !prefix.empty?
		state = DEFAULT
	end

	if !line.index("/db_xref=\"GeneID:").nil?
		id = line.strip.gsub("/db_xref=\"GeneID:","").gsub("\"","")
		mrnas[id.to_i] = true if state == MRNA
		cds[id.to_i] = true if state == CDS
		genes[id.to_i] = true if state == GENE
	end

end

genes.keys.each do |key|
	puts "#{key} : #{cds.has_key?(key)} \t: #{mrnas.has_key?(key)}"
end

puts "------------"

cds.keys.each do |key|
	puts "#{key} : #{mrnas.has_key?(key)} \t: #{genes.has_key?(key)}"
end

puts "------------"

mrnas.keys.each do |key|
	puts "#{key} : #{cds.has_key?(key)} \t: #{genes.has_key?(key)}"
end

puts "#{genes.keys.length} / #{cds.keys.length} / #{mrnas.keys.length} "