class CheckMate

require 'rubygems'
require 'csv'
require 'json'

 def json_to_csv(file_read, file_write)     
		CSV.open(file_write, 'w') do |csv|
    			JSON.parse(File.open(file_read).read).each do |k,v|
		                if v.kind_of? String
		                   csv << [k,v]
		                elsif v.kind_of? Hash
		                   csv << [k]
                                   v.each {|k,v| csv << ["",k,v]}
		                end          
			end
		end

 end

 def csv_to_json(file_read, file_write)
    
           string = IO.read(file_read)
           array = string.split("\n")
           array_splitted = array.map {|a| a.split(",")} 
           array_transpose = array_splitted.transpose
           
              h = Hash.new {|hash, key| hash[key] = []}
              for k in 0..array_transpose.size - 1 
                h[k] << array_transpose[k]
              end

   File.open(file_write, 'w') {|file| file.write(JSON.pretty_generate(h))}

 end


end

ani = CheckMate.new
ani.csv_to_json('coeffs_file.csv', 'example.json')

