set.seed(99)
files <- list.files()

n = length(files) - 1
v <- round(runif(n)*(n-1))+1

chosen.files <- files[v]

write(chosen.files,file = "chosen_files.txt")
write(v,file="chosen_files_number.txt")

