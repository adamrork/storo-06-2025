# These lists and vectors are a bit cumbersome to include in every file that requires #
# them so we will keep them in this R document and load them in where necessary #

# Contrasts for X vs. Y #
deContrasts <- c("NT_D10_vs_D0_D00",
                 "NT_D24_vs_D0_D00",
                 "PRG4_D10_vs_NT_D10",
                 "PRG4_D24_vs_NT_D24",
                 "PRG4.VEGF_D10_vs_NT_D10",
                 "PRG4.VEGF_D24_vs_NT_D24",
                 "PRG4.VEGF_D10_vs_PRG4_D10",
                 "PRG4.VEGF_D24_vs_PRG4_D24",
                 "NT_D24_vs_NT_D10",
                 "PRG4_D24_vs_PRG4_D10",
                 "PRG4.VEGF_D24_vs_PRG4.VEGF_D10")

# Create a vector of contrasts which have rather low quality PWF nullp() plots #
low.qual.pwfs <- c("PRG4_D10_vs_NT_D10_DOWN")

# Directions of differential expression #
regDirection <- c("UP", "DOWN")

# GO Categories #
goCategories <- c("BP", "MF", "CC")

# RRVGO clustering levels #
clusterLevels <- c("50", "70", "90")

# clusterProfiler methods #
annotStrategies <- c("goEnr", "kgEnr", "goGSEnr", "kgGSEnr")

# clusterProfiler source #
annotTypes <- c("GO", "KEGG")

# clusterProfiler enrichment type #
annotMethods <- c("Enr", "GSEnr")
