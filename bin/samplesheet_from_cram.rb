this_dir = Dir.pwd

root_directory = ARGV.shift

crams = Dir["#{root_directory}/*.cram"]


puts "Individual;Sample;Cram;Bai"

crams.each do |cram_file|
	file = cram_file.split("/")[-1]
	individual = file.split(".")[0]
	sample = file.split(".")[1]
	puts "#{individual};#{sample};#{File.expand_path(cram_file)};#{File.expand_path(cram_file)}.bai"
end

