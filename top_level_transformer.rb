require "date"

class TopLevelTransformer

  def initialize
  end

  def transform(data)
    result = {}
    result["LOCUS"] = parse_locus(data)
    result["ORGANISM"] = parse_organism(data)
    result["DESCRIPTION"] = data["DEFINITION"].join(" ")
    result["version"] = data["VERSION"].join(" ")
    return result
  end

  def parse_locus(data)
    locus_res = {}
    splitted = data["LOCUS"][0].split(" ")
    locus_res["refseq_id"] = splitted[0]
    locus_res["lengthh"] = splitted[1].to_i
    locus_res["date"] = parse_date(splitted[-1])
    return locus_res
  end

  def parse_organism(data)
    org_res = {}
    org_res["name"] = data["ORGANISM"][0].strip
    org_res["taxonomy"] = data["ORGANISM"][1..-1].join(";").gsub(".","").split(";")
    return org_res
  end

  def parse_date(date_str)
    splitted = date_str.split("-")
    splitted[1] = Date::ABBR_MONTHNAMES.compact.map(&:downcase).find_index(splitted[1].downcase) + 1
    return splitted.join(".")
  end

end