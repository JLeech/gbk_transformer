require_relative "top_level_transformer.rb"
require_relative "second_level_transformer.rb"

class GbkTransformer

  attr_accessor :out_file
  attr_accessor :in_file
  attr_accessor :state

  TOPLEVEL = 0
  FEATURES = 1
  ORIGIN = 2

  def initialize(gbk_path)
    File.write(gbk_path.gsub("gbk", "gbk_opt"), "")
    self.in_file = File.open(gbk_path, "r")
    self.out_file = File.open(gbk_path.gsub("gbk", "gbk_opt"), "a")
    self.state = TOPLEVEL
  end

  def transform
    top_data = Hash.new { |hash, key| hash[key] = [] }
    feat_data = []
    keys = [ "gene","source","CDS"]
    data_key = ""
    while (line = in_file.gets)
      if line.strip.start_with?("//")
        top_data = Hash.new { |hash, key| hash[key] = [] }
        data_key = ""
        self.state = TOPLEVEL
      end
      if state == TOPLEVEL
        next if line.length < 13
        prefix = line[0..11].strip
        value = line[12..-1].strip

        if prefix == "FEATURES"
          top_transformed = TopLevelTransformer.new.transform(top_data)
          top_transformed["source_file_name"] = self.in_file
          self.out_file.puts(top_transformed)
          #self.result = TopLevelTransformer.new.transform(top_data)
          self.state = FEATURES
          top_data = Hash.new { |hash, key| hash[key] = [] }
          next
        end
        data_key = prefix if !prefix.empty?
        top_data[data_key] << value
      end

      if state == FEATURES
        prefix = line[0..20].strip
        value = line[21..-1].strip if line.length > 20
        if prefix == "ORIGIN"
          self.state = ORIGIN
          next
        end
        if !prefix.empty?
          if ((keys.include?(data_key)) || data_key.end_with?("RNA") )
            self.out_file.puts(SecondLevelTransformer.new.transform(data_key, feat_data))
          end
          data_key = prefix
          feat_data = []
        end
        feat_data << value
      end
    end
  end

end



files = File.read('file_with_gbks.txt').split("\n")


files.each do |gbk_path|
  next if gbk_path.start_with?("#")
  transformer = GbkTransformer.new(gbk_path)
  transformer.transform

end

