require 'rubygems'
require 'progressbar'

def progress


in_name     = "test.gff"
out_name    = "output.fasta"

in_file     = File.new(in_name, "r")
out_file    = File.new(out_name, "w")

in_size     = File.size(in_name)
batch_bytes = ( in_size / 100 ).ceil
total       = 0
p_bar       = ProgressBar.new('Parsing', 100)

buffer      = in_file.sysread(batch_bytes)
while total < in_size do
 out_file.syswrite(buffer)
 p_bar.inc
 total += batch_bytes
 if (in_size - total) < batch_bytes
   batch_bytes = (in_size - total)
 end
 buffer = in_file.sysread(batch_bytes)
end
p_bar.finish


end

progress
