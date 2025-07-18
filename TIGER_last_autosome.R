
# Function to return the number of the last autosome based on genome build
TIGER_last_autosome <- function(genome_build) {
    num <- switch(genome_build,
                  'mm10' = 19,
                  'dm6' = 5,
                  22)  # otherwise case - default return value
    return(num)
}




# Translation Notes:
# 
# Function Definition: Translated Matlab's function num = TIGER_last_autosome(genome_build) to R's TIGER_last_autosome <- function(genome_build) format.
# 
# Switch Statement:
# 
# Matlab uses switch variable case 'value' ... otherwise syntax
# R uses switch(variable, 'value1' = result1, 'value2' = result2, default_value) syntax
# The R switch function automatically returns the matched value or the default (last unmatched value)
# Case Handling:
# 
# 'mm10' case: Returns 19 (mouse genome - 19 autosomes)
# 'dm6' case: Returns 5 (fruit fly genome - 5 autosomes)
# Default case: Returns 22 (likely human genome - 22 autosomes)
# Return Value: Added explicit return(num) for clarity, though R functions automatically return the last evaluated expression.
# 
# Variable Assignment: Used num <- to assign the switch result, maintaining the same variable name as the original Matlab code.
# 
# Functionality: The translated function maintains identical behavior - 
# it takes a genome build string and returns the corresponding number of autosomes for that species.
# 
# No dependencies or additional packages are required for this translation. The function will work identically to the original Matlab version.
# 

