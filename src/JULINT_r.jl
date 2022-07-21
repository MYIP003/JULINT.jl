#Read .int file from ISCE and StaMPS (.int format = float 32)
# XML file should be in the same folder with .int file
# Edited 12/08/2021

using LightXML
using Mmap
using LinLogQuantization
#=
function julint_r(filename::String)

   xml_filename = string(filename,".xml");
   xdoc = parse_file(xml_filename);
   xroot = root(xdoc);
   ces = collect(child_elements(xroot));   # with out name: imageFile
	
   if attribute(ces[12],"name") == "length" && attribute(ces[16],"name") == "width"
	c = find_element(ces[12],"value");
	lines = parse(Int64,content(c));
	d = find_element(ces[16],"value");
	width = parse(Int64,content(d));
   else
	println("Error, cannot find length and width from xml file.")     
   end

   data=read(filename);
   b = reinterpret(Float32, data)
   data=reshape(transpose(b),width*2,lines);
   data = transpose(data);

   INT = Array{ComplexF32}(undef,width, lines);
   INT = Complex.(data[:,1:2:end],data[:,2:2:end]);

   return INT
end
=#

function julint_r_fast(filename::String)

   xml_filename = string(filename,".xml");
   xdoc = parse_file(xml_filename);
   xroot = root(xdoc);
   ces = collect(child_elements(xroot));   # with out name: imageFile

   if attribute(ces[12],"name") == "length" && attribute(ces[16],"name") == "width"
        c = find_element(ces[12],"value");
        lines = parse(Int64,content(c));
        d = find_element(ces[16],"value");
        width = parse(Int64,content(d));
   else
        println("Error, cannot find length and width from xml file.")
   end
   
   s = open(filename);
   INT = transpose(Mmap.mmap(s, Matrix{ComplexF32},(width,lines)));
   
   ### Old script:
   #data = transpose(Mmap.mmap(s, Matrix{Float32},(width*2,lines)));
   #INT_r = INT_i = zeros(Float32,lines,width); # Pre-allocation
   #ph = A = zeros(Float32,lines,width); # Pre-alloation 
   #INT_r = data[:,[1:2:width*2;]];
   #INT_i = data[:,[2:2:width*2;]];
   #ph = atan.(INT_i,INT_r);		   # ph = angle.(INT)
   #A = sqrt.(abs2.(INT_r)+ abs2.(INT_i)); # A = abs.(INT)
   #INT = Array{ComplexF32}(undef,width, lines);
   #INT = Complex.(data[:,[1:2:width*2;]],data[:,[2:2:width*2;]]);
      
   return INT
   close(s)
end

# Read file and output phase and amplitude(isce_minrefdem.int)
# e.g. ph,A = julint_pha(filename)
# Old function
#function julint_pha(filename::String)
#
#  xml_filename = string(filename,".xml");
#   xdoc = parse_file(xml_filename);
#   xroot = root(xdoc);
#   ces = collect(child_elements(xroot));   # with out name: imageFile
#
#   if attribute(ces[12],"name") == "length" && attribute(ces[16],"name") == "width"
#        c = find_element(ces[12],"value");
#        lines = parse(Int64,content(c));
#        d = find_element(ces[16],"value");
#        width = parse(Int64,content(d));
#   else
#        println("Error, cannot find length and width from xml file.")
#   end
#   
#   s = open(filename);
#   data = transpose(Mmap.mmap(s, Matrix{Float32},(width*2,lines)));
#   
#   INT_r = data[:,[1:2:width*2;]];
#   INT_i = data[:,[2:2:width*2;]];
#   ph = atan.(INT_i,INT_r);               # ph = angle.(INT)
#   A = sqrt.(abs2.(INT_r)+ abs2.(INT_i)); # A = abs.(INT)
#    
#   return ph,A
#   close(s)
#end

## New function
function julint_pha(filename::String)

   xml_filename = string(filename,".xml");
   xdoc = parse_file(xml_filename);
   xroot = root(xdoc);
   ces = collect(child_elements(xroot));   # with out name: imageFile

   if attribute(ces[12],"name") == "length" && attribute(ces[16],"name") == "width"
        c = find_element(ces[12],"value");
        lines = parse(Int64,content(c));
        d = find_element(ces[16],"value");
        width = parse(Int64,content(d));
   else
        println("Error, cannot find length and width from xml file.")
   end
   
   s = open(filename);
   data = transpose(Mmap.mmap(s, Matrix{ComplexF32},(width,lines)));
   ph = angle.(data);
   A = abs.(data);

   return ph,A
   close(s)
end

# Function to read correlation map
function julint_r_cor(filename::String)

   xml_filename = string(filename,".xml");
   xdoc = parse_file(xml_filename);
   xroot = root(xdoc);
   ces = collect(child_elements(xroot));   # with out name: imageFile

   if attribute(ces[10],"name") == "length" && attribute(ces[14],"name") == "width"
        c = find_element(ces[10],"value");
        lines = parse(Int64,content(c));
        d = find_element(ces[14],"value");
        width = parse(Int64,content(d));
   else
        println("Error, cannot find length and width from xml file.")
   end

   s = open(filename);
   data = transpose(Mmap.mmap(s, Matrix{Float32},(width,lines)));

   return data
   close(s)
end

# Function to read merge/.int map
function julint_r_int(filename::String)
   
   xml_filename = string(filename,".xml");
   xdoc = parse_file(xml_filename);
   xroot = root(xdoc);
   ces = collect(child_elements(xroot));   # with out name: imageFile

   if attribute(ces[11],"name") == "length" && attribute(ces[15],"name") == "width"
        c = find_element(ces[11],"value");
        lines = parse(Int64,content(c));
        d = find_element(ces[15],"value");
        width = parse(Int64,content(d));
   else
        println("Error, cannot find length and width from xml file.")
   end

   s = open(filename);
   data = transpose(Mmap.mmap(s, Matrix{Float32},(width,lines)));

   return data
   close(s)
end

###########################################################################
#=
# De-quantized UInt to Float
function Base.Array{T}(n::Integer,Q::LinQuantArray) where {T<:AbstractFloat}
    Qmin = Q.min                # min as Float64
    Qmax = Q.max                # max as Float64
    Δ = (Qmax-Qmin)/(2^n-1)     # linear spacing

    A = similar(Q,T)

    @inbounds for i in eachindex(A)
        # convert Q[i]::UInt to Float64 via *
        # then to T through =
        A[i] = Qmin + Q[i]*Δ
    end

    return A
end


# De-quantized the lin8 Array to log16
Amplitude = A;
function Base.Array{T}(n::Integer,Q::LinQuantArray) where {T<:Unsigned}
    Qmin = Q.min                # min as Float64
    Qmax = Q.max                # max as Float64
    Δ = (Qmax-Qmin)/(2^n-1)     # linear spacing

    A = similar(Q,T)

    @inbounds for i in eachindex(A)
        # convert Q[i]::UInt to Float64 via *
        # then to T through =
        A[i] = Qmin + Q[i]*Δ
    end
    Amin = log(minpos(Amplitude));              # minimum of value range
    Amax = log(maximum(Amplitude));              # maximum of value range
    return LogQuantArray{T,ndims(A)}(A,Amin,Amax)
end
Base.Array{T}(Q::LinQuantArray{UInt8,N}) where {T,N} = Array{T}(8,Q)

n = 8
Q = linAlogA
Qmin = Q.min                # min as Float64
Qmax = Q.max                # max as Float64
Δ = (Qmax-Qmin)/(2^n-1)     # linear spacing
A = similar(Q,Float16)
for i in eachindex(A)
        # convert Q[i]::UInt to Float64 via *
        # then to T through =
        A[i] = Qmin + Q[i]*Δ
end
k = 1

B = similar(Q,Int64);
for i in eachindex(A)
 if A[i] == Inf16
 B[k] = i;
 k+=1;
 end
end
=#
########################################################################
