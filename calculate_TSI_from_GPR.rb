class ReadGpr

require 'rinruby'
require 'csv'

class NoGprError < StandardError
end

 #This method handles replicate .gpr files. The method reads individual replicate along with main gpr by applying
 #readGpr method. It stores the calculates TSI and names in separate arrays. The names array carries same probe names.
 #But the TSI has different intensities and the TSI array stores it as a 2D array. Each element of TSI has a calculated
 # TSI of replicate. Then the method apply transpose function to carry same probe values in a single elemnet. After this
 #step the negataive numbers are converted to zero and then airthmatic mean is applied to all elements of a individual
 #element of TSI array. if [0,0,0] the 0/3 gives us 0. if [0,0,24] then 24/3 gives 8. In other way it provides a kind of 
 # negative weight for a negative or zero value of a probe in a replicate. if [1, 21, 101] then 123/3 = 41
 def add_replicates

     #create path variable
     path = File.join("/home/jarvis","test_app","extra_scripts")

     files = Dir["#{path}/*.gpr"]
     dummynames = Array.new
     dummyTSI   = Array.new

     files.each_with_index do |file, i|
       dummynames[i], dummyTSI[i] = readGpr(file)
     end   

     tsi_container = dummyTSI.transpose
     tsi_container.map! {|e| e.map! {|f| f < 0? 0 : f}}
     tsi_container.map! {|e| e.inject(:+)/e.size}

     names_container = dummynames[0]

       
        
    puts tsi_container.to_s
    puts names_container.to_s

    puts tsi_container.size
    puts names_container.size

 end 

 def readGpr(file_path)
  
  read_array, header_removed = [], []
   begin
	     read = IO.binread(file_path)
             read = read.encode('UTF-8', :invalid => :replace, :undef => :replace)
	    
		     if read.valid_encoding?
			 read_array = read.split("\n")         
		     else
			 read_array = read.encode!("ASCII-8BIT","ASCII-8BIT", invalid: :replace).split("\n")
		     end
	    
	      mod_array = read_array.map {|e| e.split("\t")}  
	      
	      element_stabilized = mod_array.map {|element| element.join(",").gsub("\"","").split(",")} 

		      if element_stabilized[0].include?("ATF")
			 header_removed = element_stabilized.drop_while {|i| i unless i.include?("Block")}
		       else
			 raise NoGprError, "File does not seem to be gpr formatted. Check the file"
		      end

              column_based_array = header_removed.transpose
              @name, @dia, @f633_mean, @b633_mean = getColumns(column_based_array)
              @probeNames, @tsiList = calTotalSignalIntensity(@name, @dia, @f633_mean, @b633_mean)              
              #@probeNames, @tsiList = sortGprTsiList(@name, @get_tsi_list)
              #puts "#{@probeNames}"
              #puts "#{@tsiList}"

              #dummy = Array.new
              #for i in 0..@probeNames.size              
              #   dummy.push("#{@probeNames[i]}\t#{@tsiList[i]}")
              #end
              
              #dummy2 = dummy.map {|e| e.split.to_a}           
              #send_file("output_tsi_list.csv", dummy2)

              return @probeNames, @tsiList
               
    rescue Exception => e
              puts "something went wrong :##" + e.message.to_s + e.backtrace.to_s
    end 

 end 

 def getColumns(array=[])
     name, dia, f633_mean, b633_mean = [], [],[],[]
	   begin
		     array.map do |element|     
			       case
				       when element.include?("Name") then name << element
				       when element.include?("Dia.") then dia << element
				       when element.include?("F633 Mean") then f633_mean << element
				       when element.include?("B633 Mean") then b633_mean << element       
			       end
		     end
	   
	    rescue Exception => e
		   e.message
                   e.backtrace.inspect
	    end

    return name, dia, f633_mean, b633_mean 
 end

	 def calTotalSignalIntensity(probeNameList, diameter, forground, background)

	  names = probeNameList.flatten
	    names.shift
	    filterNames = names.uniq 
	    
	    dia = diameter.flatten
	    dia.shift
	    f633 = forground.flatten
	    f633.shift
	    b633 = background.flatten
	    b633.shift   
	     
	     names = partition_array(names)
	     counts = names.length

          #puts "#{names[0]}"

	   R.assign "counts", counts
	   for i in 1..names.count
	      R.assign "name#{i}", names[i-1]
	   end

     #Formula for calculating Total Signal Intensity
     #(F633_mean - B633_mean)*3.14*diameter^2*1/4
     dia = partition_array(dia)
     for i in 1..dia.count
      R.assign "dia#{i}", dia[i-1]
     end
     #R.assign "dia", dia

     f633 = partition_array(f633)
     for i in 1..f633.count
      R.assign "f633#{i}", f633[i-1]
     end
     #R.assign "f633", f633
     
     b633 = partition_array(b633)
     for i in 1..b633.count
      R.assign "b633#{i}", b633[i-1]
     end
     #R.assign "b633", b633

	 
	 R.eval <<-EOF

	 mergeVectors <- function(array, counts) {
      for (i in c(1:counts)) {
         if (i == 1) { dummy <- c(get(paste0(array,i))) } 
         else { dummy <- c(dummy, get(paste0(array,i))) }    
       }
      return(dummy)
    }
    
  names <- mergeVectors("name", counts)
  names <- as.character(names)
  dia <- mergeVectors("dia", counts)
  f633 <- mergeVectors("f633", counts)
  b633 <- mergeVectors("b633", counts)

       calTSI <- function(dia, f633, b633) {

	  dia <- as.numeric(dia)
	  f633 <- as.numeric(f633)
	  b633 <- as.numeric(b633)

	  tsi <- (f633 - b633) * 3.14 * dia * dia * 1/4

	  return(tsi)
	} 

       totalSignalIntensities <- calTSI(dia, f633, b633)

       names <- as.vector(names)
       totalSignalIntensities <- as.numeric(totalSignalIntensities)

	tab <- cbind(Name=names, F633=totalSignalIntensities)
	tab <- data.frame(tab)	
 
	allProbes <- as.character(tab[,1])
	uniqueProbeVec <- unique(allProbes) 
        uniqueProbeVecFilter <- gsub("\357\277\275\357\277\275\357\277\275M", "", uniqueProbeVec)
        print(uniqueProbeVecFilter) 

	meanTSI <- list()
	myData <- list()

	for (i in c(1:length(uniqueProbeVec))) {
	    
		myData[[i]] <- subset(tab, uniqueProbeVec[i] == tab[ , 1])
	} 

	for (j in c(1:length(uniqueProbeVec))) {

                newVec <- as.numeric(as.character(myData[[j]][, 2]))
                replicate <- as.numeric(length(newVec))

		meanTSI[[j]] <- sum(newVec)/replicate
                

	}

	meanTSI <- unlist(meanTSI)

        
      EOF
	 
	   #passing non UTF-8 char from R to ruby and vice versa throws an error... 
	    #"Error in nchar(var) : invalid multibyte string 1". 
	    #here is a workaround.
	    #Sys.setlocale('LC_ALL','C')
	    #http://stackoverflow.com/questions/6669911/error-in-nchar-when-reading-in-stata-file-in-r-on-mac

	   #Splitting with in R
	   #split(a, ceiling(seq_along(a)/3))

          
            #tsi = R.pull("totalSignalIntensities")
           tsiList = R.pull("meanTSI")
           return filterNames, tsiList


	 end


	 def partition_array(array=[], size=500)
		dummy = []
	     begin
		if array.empty?
		   raise "Input array is empty!!"
		else
		array.each_slice(size) {|element| dummy.push(element)}
		end
	     rescue Exception => e
		   e.message
	     end     

	   return dummy 
	 end

   def send_file(filename, arg=[])
     #data = args.join(',').split(',')
     file = CSV.generate do |line|
        arg.each do |element|
        line << element
        end
     end

   File.open(filename, "w+") {|e| e.puts(file)}
  end


end

readgpr = ReadGpr.new
#readgpr.readGpr("P_Chienti_4_10061160.gpr")
readgpr.add_replicates
