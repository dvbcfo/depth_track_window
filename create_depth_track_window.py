#!usr/bin/env python

import sys, string

if (len(sys.argv) == 1):
   print "Run as:\n"
   print "python create_depth_track_window.py ifolder ifile ofolder win_size [-zig] [-hg18]\n"
   print "- ifolder:  folder containing the input file"
   print "- ifile:    the input file (a coordinate-sorted bed file)"
   print "- ofolder:  folder for output files"
   print "- win_size: size of the sliding window"
   print "- [-zig]  : optional parameter, if \"-zig\" is present as an argument, bases with zero coverage will be ignored during the computation of mean genomic sequencing depth"
   print "- [-hg18]:  optional parameter, if \"-hg18\" is present as an argument, the script works with chromosome sizes derived from human genomic build hg18\n"
   sys.exit(0)

# input folder - folder containing the input file
ifolder = sys.argv[1]

if (ifolder[-1:] != "/"):
   ifolder = ifolder + "/"

# input file - a coordinate-sorted bed file
ifile = open(ifolder + sys.argv[2], "r")

# output folder - folder in which the output files should be stored
ofolder = sys.argv[3]

if (ofolder[-1:] != "/"):
   ofolder = ofolder + "/"

# sliding window size
win_size = int(sys.argv[4])

chrom_sizes = {}

# hg19 chromosome sizes
chrom_sizes["chr1"] = 249250621
chrom_sizes["chr2"] = 243199373
chrom_sizes["chr3"] = 198022430
chrom_sizes["chr4"] = 191154276
chrom_sizes["chr5"] = 180915260
chrom_sizes["chr6"] = 171115067
chrom_sizes["chr7"] = 159138663
chrom_sizes["chr8"] = 146364022
chrom_sizes["chr9"] = 141213431
chrom_sizes["chr10"] = 135534747
chrom_sizes["chr11"] = 135006516
chrom_sizes["chr12"] = 133851895
chrom_sizes["chr13"] = 115169878
chrom_sizes["chr14"] = 107349540
chrom_sizes["chr15"] = 102531392
chrom_sizes["chr16"] = 90354753
chrom_sizes["chr17"] = 81195210
chrom_sizes["chr18"] = 78077248
chrom_sizes["chr19"] = 59128983
chrom_sizes["chr20"] = 63025520
chrom_sizes["chr21"] = 48129895
chrom_sizes["chr22"] = 51304566
chrom_sizes["chrX"] = 155270560
chrom_sizes["chrY"] = 59373566
chrom_sizes["chrM"] = 16571

if ("-hg18" in sys.argv):
   chrom_sizes["chr1"] = 247249719
   chrom_sizes["chr2"] = 242951149
   chrom_sizes["chr3"] = 199501827
   chrom_sizes["chr4"] = 191273063
   chrom_sizes["chr5"] = 180857866
   chrom_sizes["chr6"] = 170899992
   chrom_sizes["chr7"] = 158821424
   chrom_sizes["chr8"] = 146274826
   chrom_sizes["chr9"] = 140273252
   chrom_sizes["chr10"] = 135374737
   chrom_sizes["chr11"] = 134452384
   chrom_sizes["chr12"] = 132349534
   chrom_sizes["chr13"] = 114142980
   chrom_sizes["chr14"] = 106368585
   chrom_sizes["chr15"] = 100338915
   chrom_sizes["chr16"] = 88827254
   chrom_sizes["chr17"] = 78774742
   chrom_sizes["chr18"] = 76117153
   chrom_sizes["chr19"] = 63811651
   chrom_sizes["chr20"] = 62435964
   chrom_sizes["chr21"] = 46944323
   chrom_sizes["chr22"] = 49691432
   chrom_sizes["chrX"] = 154913754
   chrom_sizes["chrY"] = 57772954
   chrom_sizes["chrM"] = 16571

   print "Optional parameter \"-hg18\" recognized: the script will be working with chromosome sizes derived from human genomic build hg18.\n"
else:
   print "Optional parameter \"-hg18\" not recognized: the script will be working with chromosome sizes derived from human genomic build hg19.\n"

pos_vals = {}   # relevant positions on currently processed chromosome and their seq. depth
final_pos_vals = {}   # values for positions that can no longer be updated (ready for ouput if all values in given sliding window are available)
out_FPV = 0   # last output position (passed from final_pos_vals to the output)
out_sum_FPV = 0   # sum of sequencing depths for all bases inside the sliding window
elem_count_FPV = 0   # number of bases inside the sliding window
max_FPV = 0   # current maximum position in the final_pos_vals directory

val_distrib = {}   # abundances of sequencing depth values

last_chrom = ""   # chromosome of the last processed read
last_pos = 0   # last genomic position passed from pos_vals to final_pos_vals

# an output buffer and its maximal allowed size
out_buffer = ""
out_buffer_size = 0
out_buffer_max_size = 1024*1024*16

ofile = ""

zero_pos = 0
non_zero_pos = 0
cov_sum = 0

# process the individual BED file records
for line in ifile:
   line = string.strip(line)
   line_s = line.split()

   # extract the record chromosome, the start- and the end- positions
   chrom = line_s[0]
   s_pos = int(line_s[1]) + 1
   e_pos = int(line_s[2])

   # upon chromosome change:
   # empty pos_vals and process any trailing zero-depth positions, change the output file
   if not (chrom == last_chrom):

      # skip this block if the first record of the BED file is being processed
      if not (last_chrom == ""):

         if (len(out_buffer) > 0):
            ofile.write(out_buffer)
            out_buffer = ""
            out_buffer_size = 0

         # process all positions not yet included into the final_pos_vals directory
         for i in range(last_pos + 1, chrom_sizes[last_chrom] + 1):

            # for positions with non-zero coverage
            if (pos_vals.has_key(i)):
               final_pos_vals[i] = pos_vals[i]
               max_FPV = i

               # update the value distribution directory
               if not (val_distrib.has_key(pos_vals[i])):
                  val_distrib[pos_vals[i]] = 1
               else:
                  val_distrib[pos_vals[i]] += 1

               non_zero_pos += 1
               cov_sum += pos_vals[i]

               del(pos_vals[i])

            # for positions with zero coverage
            else:
               max_FPV = i

               # update the value distribution directory
               if not (val_distrib.has_key(0)):
                  val_distrib[0] = 1
               else:
                  val_distrib[0] += 1

               zero_pos += 1         

         # process all positions not yet included in the output
         for i in range(out_FPV + 1, chrom_sizes[last_chrom] + 1):

            # at the beginning of a chromosome
            # - the effective sliding window size equals half of its potential maximum
            # - the sum of coverages within the sliding window is based on positions 1 - (win_size/2) 
            if (i == 1):
               for j in range(1, win_size/2 + 1):
                  elem_count_FPV += 1
                  if (final_pos_vals.has_key(j)):
                     out_sum_FPV += final_pos_vals[j]

            # at chromosomal positions where the first base of the sliding window changes (i.e., at all sites at least win_size/2 bases away from the beginning of the chromosome)
            # - adjust the current coverage sum value and the current window size by disregarding the newly excluded position
            if (i - win_size/2 > 1):
               elem_count_FPV -= 1

               if (final_pos_vals.has_key(i - win_size/2 - 1)):
                  out_sum_FPV -= final_pos_vals[i - win_size/2 - 1]
                  del(final_pos_vals[i - win_size/2 - 1])

            # at chromosomal positions where the last base of the sliding window changes (i.e., at all sites at least win_size/2 bases away from the end of the chromosome)
            # - adjust the current coverage sum value and the current window size by taking into account the newly included position
            if (i + win_size/2 <= chrom_sizes[last_chrom]):
               elem_count_FPV += 1

               if (final_pos_vals.has_key(i + win_size/2)):
                  out_sum_FPV += final_pos_vals[i + win_size/2]

            outval = round(float(out_sum_FPV)/elem_count_FPV, 4)

            # values used as ploidy values for iffBuilder need to be limited to < 65.535
            # due to subsequent normalization and iffBuilder's ability to take care of the rounding, this part of the code is now obsolete
            #if (outval > 65.53):
            #   outval = 65.53
            ofile.write(str(outval) + "\n")

         # re-initiate all values before processing BED records of the next chromosome
         final_pos_vals = {}
         out_FPV = 0
         out_sum_FPV = 0
         elem_count_FPV = 0
         max_FPV = 0

         pos_vals = {}
         last_pos = 0

      if not (ofile == ""):
         ofile.close()
      ofile = open(ofolder + sys.argv[2] + "." + chrom, "w")
      print("Creating a smooth coverage profile for chromosome " + chrom + " ..\n")
      last_chrom = chrom

   # output values for genommic positions that can't get updated anymore
   if (s_pos > last_pos + 1):

      # process all positions that cannot get more coverage support from input BED file's remaining records
      for i in range(last_pos + 1, s_pos):


         # for positions with non-zero coverage
         if (pos_vals.has_key(i)):
            final_pos_vals[i] = pos_vals[i]
            max_FPV = i

            # update the value distribution directory
            if not (val_distrib.has_key(pos_vals[i])):
               val_distrib[pos_vals[i]] = 1
            else:
               val_distrib[pos_vals[i]]	+= 1

            non_zero_pos += 1
            cov_sum += pos_vals[i]

            del(pos_vals[i])

         # for positions with zero coverage
         else:
            max_FPV = i

            # update the value distribution directory
            if not (val_distrib.has_key(0)):
               val_distrib[0] = 1
            else:
               val_distrib[0] += 1

            zero_pos += 1

      last_pos = s_pos - 1

      # output all positions for which none of the bases within the sliding window can get more coverage support from input BED file's remaining records
      if (out_FPV + 1 + win_size/2 <= max_FPV):
         for i in range(out_FPV + 1, max_FPV - win_size/2 + 1):

            # at the beginning of a chromosome
            # - the effective sliding window size equals half of its potential maximum
            # - the sum of coverages within the sliding window is based on positions 1 - (win_size/2) 
            if (i == 1):
               for j in range(1, win_size/2 + 1):
                  elem_count_FPV += 1
                  if (final_pos_vals.has_key(j)):
                     out_sum_FPV += final_pos_vals[j]

            # at chromosomal positions where the first base of the sliding window changes (i.e., at all sites at least win_size/2 bases away from the beginning of the chromosome)
            # - adjust the current coverage sum value and the current window size by disregarding the newly excluded position
            if (i - win_size/2 > 1):
               elem_count_FPV -= 1

               if (final_pos_vals.has_key(i - win_size/2 - 1)):
                  out_sum_FPV -= final_pos_vals[i - win_size/2 - 1]
                  del(final_pos_vals[i - win_size/2 - 1])

            # at chromosomal positions where the last base of the sliding window changes (i.e., at all sites at least win_size/2 bases away from the end of the chromosome)
            # - adjust the current coverage sum value and the current window size by taking into account the newly included position
            if (i + win_size/2 <= chrom_sizes[last_chrom]):
               elem_count_FPV += 1

               if (final_pos_vals.has_key(i + win_size/2)):
                  out_sum_FPV += final_pos_vals[i + win_size/2]

            outval = round(float(out_sum_FPV)/elem_count_FPV, 4)

            # values used as ploidy values for iffBuilder need to be limited to < 65.535
            # due to subsequent normalization and iffBuilder's ability to take care of the rounding, this part of the code is now obsolete
            #if (outval > 65.53):
            #   outval = 65.53
            out_buffer += str(outval) + "\n"
            out_buffer_size += 1

            if (out_buffer_size > out_buffer_max_size):
               ofile.write(out_buffer)
               out_buffer = ""
               out_buffer_size = 0

         out_FPV = max_FPV - win_size/2

   # update the pos_vals dictionary (increment for genomic positions s_pos to e_pos)
   for i in range(s_pos, e_pos + 1):
      if not (pos_vals.has_key(i)):
         pos_vals[i] = 1
      else:
         pos_vals[i] += 1

ifile.close()

# empty pos_vals and process trailing zero-depth positions for the last chromosome
if (len(out_buffer) > 0):
   ofile.write(out_buffer)
   out_buffer = ""
   out_buffer_size = 0

# process all positions not yet included into the final_pos_vals directory
for i in range(last_pos + 1, chrom_sizes[last_chrom] + 1):

   # for positions with non-zero coverage
   if (pos_vals.has_key(i)):
      final_pos_vals[i] = pos_vals[i]
      max_FPV = i

      # update the value distribution directory
      if not (val_distrib.has_key(pos_vals[i])):
         val_distrib[pos_vals[i]] = 1
      else:
         val_distrib[pos_vals[i]] += 1

      non_zero_pos += 1
      cov_sum += pos_vals[i]

      del(pos_vals[i])

   # for positions with zero coverage
   else:
      max_FPV = i

      # update the value distribution directory
      if not (val_distrib.has_key(0)):
         val_distrib[0] = 1
      else:
         val_distrib[0] += 1

      zero_pos += 1         

# process all positions not yet included in the output
for i in range(out_FPV + 1, chrom_sizes[last_chrom] + 1):

   # at the beginning of a chromosome
   # - the effective sliding window size equals half of its potential maximum
   # - the sum of coverages within the sliding window is based on positions 1 - (win_size/2)
   if (i == 1):
      for j in range(1, win_size/2 + 1):
         elem_count_FPV += 1
         if (final_pos_vals.has_key(j)):
            out_sum_FPV += final_pos_vals[j]

   # at chromosomal positions where the first base of the sliding window changes (i.e., at all sites at least win_size/2 bases away from the beginning of the chromosome)
   # - adjust the current coverage sum value and the current window size by disregarding the newly excluded positiona
   if (i - win_size/2 > 1):
      elem_count_FPV -= 1

      if (final_pos_vals.has_key(i - win_size/2 - 1)):
         out_sum_FPV -= final_pos_vals[i - win_size/2 - 1]
         del(final_pos_vals[i - win_size/2 - 1])

   # at chromosomal positions where the last base of the sliding window changes (i.e., at all sites at least win_size/2 bases away from the end of the chromosome)
   # - adjust the current coverage sum value and the current window size by taking into account the newly included position
   if (i + win_size/2 <= chrom_sizes[last_chrom]):
      elem_count_FPV += 1

      if (final_pos_vals.has_key(i + win_size/2)):
         out_sum_FPV += final_pos_vals[i + win_size/2]

   outval = round(float(out_sum_FPV)/elem_count_FPV, 4)

   # values used as ploidy values for iffBuilder need to be limited to < 65.535
   # due to subsequent normalization and iffBuilder's ability to take care of the rounding, this part of the code is now obsolete
   #if (outval > 65.53):
   #   outval = 65.53
   ofile.write(str(outval) + "\n")

final_pos_vals = {}
out_FPV = 0
out_sum_FPV = 0
elem_count_FPV = 0
max_FPV = 0

last_pos = 0

ofile.close()

# output the depth distribution values
dist_values = val_distrib.keys()
dist_values.sort()

genome_length = 0

for chrom in chrom_sizes:
   genome_length += chrom_sizes[chrom]

dfile = open(ofolder + sys.argv[2] + ".depth_value_distribution", "w")
dfile.write("# assembly: hg19\n")
dfile.write("# total genome length: " + str(genome_length) + " bp\n")
dfile.write("# column seq_depth: sequencing depth\n")
dfile.write("# column bp_count: number of genomic locations (1-bp sites) with respective sequencing depth\n")
dfile.write("# column genome_fraction_%: column bp_count expressed as fraction of the genome length\n")
dfile.write("# column acc_genome_fraction_%: accumulation of values from column genome_fraction_%\n")
dfile.write("seq_depth\tbp_count\tgenome_fraction_%\tacc_genome_fraction_%\n")

acc_b_count = 0

for value in dist_values:
   b_count = val_distrib[value]
   acc_b_count += b_count
   g_fraction = round(float(b_count)*100/genome_length, 5)
   acc_g_fraction = round(float(acc_b_count)*100/genome_length, 5)

   dfile.write(str(value) + "\t" + str(b_count) + "\t" + str(g_fraction) + "\t" + str(acc_g_fraction) + "\n")

dfile.close()

print("Number of zero-coverage positions: " + str(zero_pos))
print("Number of non-zero-coverage positions: " + str(non_zero_pos))
print("Total number of genomic positions: " + str(zero_pos + non_zero_pos) + " (check: " + str(genome_length) + ")")
print("Total coverage sum: " + str(cov_sum))
print("Genomic mean (zero-coverage positions considered): " + str(round(float(cov_sum)/(zero_pos+non_zero_pos), 4)))
print("Genomic mean (zero-coverage positions not considered): " + str(round(float(cov_sum)/non_zero_pos, 4)) + "\n")

gen_mean = float(cov_sum)/(zero_pos+non_zero_pos)

if ("-zig" in sys.argv):
   gen_mean = float(cov_sum)/non_zero_pos
   print "Optional parameter \"-zig\" recognized: bases with zero coverage are not considered in the genomic mean computation.\n"
else:
   print "Optional parameter \"-zig\" not recognized: bases with zero coverage are considered in the genomic mean computation.\n"

# creation of chromosome-wise wiggle-files based on the window-based BED read support and genomic coverage mean
for chrom_name in chrom_sizes:
   print("Creating a wiggle file for chromosome " + chrom_name + " ..\n")
   cdfile = open(ofolder + sys.argv[2] + "." + chrom_name, "r")
   wigfile = open(ofolder + sys.argv[2] + "." + chrom_name + ".wig", "w")

   out = False
   start_offset = 1

   for inline in cdfile:

      inline = string.strip(inline)

      if ((not out) and (float(inline) == 0.0)):
         start_offset += 1

      elif ((not out) and (float(inline) != 0.0)):
         out = True
         wigfile.write("fixedStep chrom=" + chrom_name + " start=" + str(start_offset) + " step=1\n")
         wigfile.write(str(round(float(inline)/gen_mean, 4)) + "\n")

      else:
         wigfile.write(str(round(float(inline)/gen_mean, 4)) + "\n")

   cdfile.close()
   wigfile.close()

print "All done.\n"
