##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  any later version.                                                    #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program. If not, see <http://www.gnu.or/licenses/>.   #
#                                                                        #
##########################################################################

## consistent error / warning messages;
..msg <- list(error =
              c(non.dna = "Function only accepts ShortReadQ or DNAStringSet objects",
                empty = "Data object is empty",
                no_file = "Output filename is missing"
                ), warn = c()
              )

phamming.distance <- function (data, output_filename) {

  objectType <- class(data)
  if(!length(data)) stop(..msg$error["empty"])

  if (is.null(output_filename)) stop(..msg$error["empty"])
  
  if (objectType=='ShortReadQ') {  
    data <- sread(data)
  } else if (objectType!='DNAStringSet') {
    stop(..msg$error["non.dna"])
  }

  sample_width <- width(data[1])
  number_of_samples <- length(data)

  if(sample_width<1 || number_of_samples<2) stop(..msg$error["empty"])

  return_val <- .C("phamming",
                   as.character(IRanges::unlist(data)),
                   as.character(output_filename),
                   as.integer(sample_width),
                   as.integer(number_of_samples)                   
                   )

  # Return values from the interface have meaning.
  #  0    -->     success
  # -1    -->     MPI is not initialized
  
  return(return_val)
}
