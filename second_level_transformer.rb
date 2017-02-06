class SecondLevelTransformer

  attr_accessor :result

  def initialize
    self.result = {}
  end

  def transform(data_key, data)
    parse_gene(data) if data_key == "gene"
    parse_source(data) if data_key == "source"
    parse_cds(data) if data_key == "CDS"
    if data_key == "mRNA"
      parse_mrna(data) 
    elsif data_key.end_with?("RNA")
      parse_rna(data)
    end
    return result
  end

  def parse_gene(data)
    self.result["start"], self.result["end"], self.result["backward_chain"],_,_  = parse_range(data)
    attributes = parse_attributes(data)
    self.result["name"] = attributes["gene"]
    self.result["is_presudo"] = !(attributes.keys & ["\"pseudo\"","\"pseudogene\""]).empty?
    self.result["type"] = "gene"
  end

  def parse_source(data)
    attributes = parse_attributes(data)
    self.result["db_mitochondria"] = (attributes["organelle"] == "\"mitochondrion\"")
    self.result["taxonomy_xref"] = attributes["db_xref"]
    self.result["type"] = "source"
  end

  def parse_cds(data)
    st, fin, complement, starts, fins  = parse_range(data)
    self.result["isoform"] = {"cdsStart" => st,"cdsEnd" => fin,"exonsCdsCount" => starts.length}
    self.result["gene"] = {"isProteinButNotRna" => true, "hasCDS" => true, "startCode" => st, "endCode" => fin}
    self.result["type"] = "CDS"

    attributes = parse_attributes(data)

    ["protein_id","db_xref","product","note"].each do |cur|
      self.result["isoform"][cur] = attributes[cur] if attributes.has_key?(cur)
    end
    self.result["exons"], self.result["introns"] = make_exons_introns(complement, starts, fins)
  end

  def parse_mrna(data)
    self.result["mrna_start"], self.result["mrna_end"], _, self.result["starts"], self.result["fins"]  = parse_range(data)
    attributes = parse_attributes(data)
    self.result["type"] = "MRNA"
    self.result["exons_mrna_count"] = self.result["starts"].length
    self.result["mrna_length"] = self.result["mrna_end"] - self.result["mrna_start"] + 1
    self.result["protein_id"] = attributes["protein_id"]
    self.result["protein_xref"] = attributes["db_xref"]
    self.result["product"] = attributes["product"]
    self.result["note"] = attributes["note"]
  end

  def parse_rna(data)
    self.result["is_rna"] = true
    self.result["type"] = "RNA"
  end

private

  def parse_range(data)
    range_strs = ""
    data.each do |part|
      break if part[0] == "/"
      range_strs += part.gsub("\n","")
    end
      
    complement = range_strs.start_with?("complement")

    range_strs.gsub!("complement(","")
    range_strs.gsub!("join(","")
    range_strs.gsub!(/(<|\(|\)|>)/,"")

    puts range_strs.split(",")

    starts = []
    fins = []
    range_strs.split(",").each do |rang| 
      splt = rang.split("..")
      starts << splt[0].to_i
      fins << splt[1].to_i
    end

    return [ starts.min, fins.max, complement, starts, fins]
  end

  def parse_attributes(data)

    attributes = Hash.new { |hash, key| hash[key] = "" }
    on_place = false

    cur_attr = ""
    data.each do |part|
      on_place = true if part[0] == "/"
      if on_place
        if part[0] == "/"
          if part.index("=").nil?
            attributes[part[1..-1]] << ""
          else
            splitted = part.split("=")
            cur_attr = splitted[0][1..-1]
            attributes[cur_attr] << splitted[1]
          end
        else
          attributes[cur_attr] << part
        end
      end
    end
    return attributes
  end

  def make_exons_introns(complement, starts, fins)

    return [] if starts.length == 0
    
    exons = []
    introns = []

    start_index = complement ? (starts.length - 1) : 0
    fin_index = complement ? -1 : starts.length 
    inc = complement ? -1 : 1

    phase = 0

    index = start_index
    while(index != fin_index)
      start = starts[index]
      fin = fins[index]
      exon = {}
      exon["start"] = start
      exon["end"] = fin
      exon["start_phase"] = phase
      exon["end_phase"] = (phase + fin - start + 1)%3
      phase = exon["end_phase"]
      exons << exon
      index += inc
    end

    if(exons.length == 1)
      exon = exons.first
      exon["index"] = 0
      exon["rev_index"] = 0
      exon["type"] = "one_exon"
    else
      exons.each_with_index do |exon, index|
        exon["index"] = index
        exon["rev_index"] = exons.length - index - 1
        if (index == 0)
          exon["type"] = "start" 
        elsif ( exons.length-1 == index )
          exon["type"] = "end"
        else
          exon["type"] = "inner"
        end
        if index > 0
          intron = {}
          prev_exon = exons[index-1]
          intron["start"] = complement ? exon["end"] + 1 : prev_exon["end"]+1
          intron["end"] = complement ? prev_exon["start"] - 1 : prev_exon["start"]-1
          intron["index"] = index - 1
          intron["rev_index"] = exons.length - index - 1
          intron["phase"] = prev_exon["end_phase"]
          intron["length_phase"] = (intron["end"] - intron["start"] + 1)%3
          introns << intron
        end
      end
    end

    return [exons, introns]

  end

end