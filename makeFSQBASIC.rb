#!/usr/bin/ruby
#---

require 'getoptlong'



now = Time.now                    

puts now



# ///////////////////////////////////////////////////////////////

# ///////////////////////////////////////////////////////////////


 result1 = %x{gcc -lrt -lasound -lpthread -lm -lwiringPi -o fsqBasic fsqBasic.c dspmath.c libportaudio.a}



 puts "<< build complete >>"
#  result1 = %x{mono file1.exe #{ARGV[0].to_i}}  






  

