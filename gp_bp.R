# Load necessary libraries
## Osei Asibey - oo324@exeter.ac.uk
## 2024.12.18 -- Load UK Biobank subset of gp_clinical that match some search terms
##               assumes exported file is in RAP project directory /ukbrapr_data/
##               It is worth double-checking the matched rows look "right" -- have accidental extra codes been caught by the grep?

#Setting Directory                      
setwd("/mnt/project/projects/Osei/")

# install required packages if missing
install.packages(setdiff(c("readr","stringr"), rownames(installed.packages())))
library(stringr)
library(dplyr)



# Step 1: Read the table into R
# Assuming the table is saved as a CSV file named "terms_table.csv"
terms_table <- readr::read_csv("/mnt/project/projects/Osei/readv2_v3.20241203.cvs")

# Extract the relevant columns (e.g., "READV2_CODE" and "TERMV3_CODE")
search_terms <- unique(c(terms_table$READV2_CODE, terms_table$TERMV3_CODE))

# Remove NA values (if any)
search_terms <- search_terms[!is.na(search_terms)]

# Combine the terms into a single string separated by '|'
flattened_terms <- str_flatten(search_terms, collapse = "|")

# Construct the grep command
search_string <- str_c("grep -E ", sprintf('"%s"', flattened_terms), " /mnt/project/ukbrapr_data/gp_clinical.tsv")

# Output the search string
cat(search_string)

# get variables names of file
headers <- colnames(readr::read_tsv("/mnt/project/ukbrapr_data/gp_clinical.tsv", n_max=1, show_col_types=FALSE, progress=FALSE))

# load subset of gp_scripts rows that match one of the search terms above
gp_bp <- readr::read_tsv(pipe(search_string), col_names=headers)


# Create BP Categories
gp_bp <- gp_bp |>
  mutate(
    bp_cat1 = case_when(
      read_2 %in% c("2469.", "246b.", "246d.",  "246e.", "246l.", "246N.", "246n1", "246o0", "246Q.", "246S.", "246W.", "246Y.") ~ 1, 
      read_2 %in% c("246a.", "246A.", "246c.", "246f.", "246m.", "246n0", "246o1", "246P.", "246R.", "246T.", "246V.", "246X.") ~ 2,
      read_2 %in% c("2469", "246..", "246C.", "246D.", "246E.", "246g.", "246h.", "246H.", "246i.", "246j.", "246K.", "246L.", "246n.", "246o.") ~ 3,
      TRUE ~ NA_real_),
    bp_cat2 = case_when(
      read_3 %in% c("2469.", "246b.", "246d.",  "246e.", "246l.", "246N.", "246n1", "246o0", "246Q.", "246S.", "246W.", "246Y.") ~ 1, 
      read_3 %in% c("246a.", "246A.", "246c.", "246f.", "246m.", "246n0", "246o1", "246P.", "246R.", "246T.", "246V.", "246X.") ~ 2,
      read_3 %in% c("2469", "246..", "246C.", "246D.", "246E.", "246g.", "246h.", "246H.", "246i.", "246j.", "246K.", "246L.", "246n.", "246o.") ~ 3,
      TRUE ~ NA_real_)
  )

# Combining gp_cat from readv2 and readv3
gp_bp$bp_cat <- gp_bp$bp_cat1
gp_bp$bp_cat <- ifelse(is.na(gp_bp$bp_cat), gp_bp$bp_cat2, gp_bp$bp_cat) 



# Define the file name and target path
filename <- "gp_bp"
target_path <- "/mnt/project/projects/Osei/gp_bp.R"

# Construct the dx upload command
dx_command <- sprintf("dx upload %s --path %s", filename, target_path)

# Execute the command
system(dx_command)
