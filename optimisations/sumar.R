setwd(dir = "./array_csvs/")
array_fles = sort(list.files())
setwd("../")

setwd(dir = "./main_outputs/")
main_fles = sort(list.files())
setwd("../")

setwd(dir = "./log_files/")
log_fles = sort(list.files())
setwd("../")


data = data.frame(name = array_fles,
                  arrays.no = rep(0, length(array_fles)),
                  array.tot.size = rep(0, length(array_fles)),
                  arrays.with.repres.no = rep(0, length(array_fles)),
                  array.with.repres.tot.size = rep(0, length(array_fles)),
                  arrays.tot.representative.size = rep(0, length(array_fles)),
                  repeat.no = rep(0, length(array_fles)),
                  repeat.tot.len = rep(0, length(array_fles)),
                  out_6 = rep(0, length(array_fles)), #6 to 7
                  out_7 = rep(0, length(array_fles)), #7 to 8
                  out_8 = rep(0, length(array_fles)), 
                  out_9 = rep(0, length(array_fles)), #9 to 11
                  out_11 = rep(0, length(array_fles)),
                  f8_t0 = rep(0, length(array_fles)),
                  f8_t1 = rep(0, length(array_fles)),
                  f8_tk = rep(0, length(array_fles)))

for(i in 1 : length(array_fles))
{
  print(i)
  arrays = read.csv(file = file.path("array_csvs", array_fles[i]))
  arrays$width = arrays$end - arrays$start + 1
  data$arrays.no[i] = nrow(arrays)
  data$array.tot.size[i] = sum(arrays$width)
  data$arrays.with.repres.no[i] = nrow(arrays[arrays$representative != "", ])
  data$array.with.repres.tot.size[i] = sum(arrays$width[arrays$representative != ""])
  data$arrays.tot.representative.size[i] = sum(nchar(arrays$representative))
  data$repeat.no[i] = sum(arrays$repeats_number)
  data$repeat.tot.len[i] = sum(arrays$repeats_number * arrays$median_repeat_width)

  main_log = readLines(file.path("main_outputs", main_fles[i]))
  for(j in 1 : length(main_log)) {
    if(grepl("06 / 13", main_log[j])) {
      data$out_6[i] = as.numeric(strsplit(main_log[j], split = " ")[[1]][length(strsplit(main_log[j], split = " ")[[1]])])
    }
    if(grepl("07 / 13", main_log[j])) {
      data$out_7[i] = as.numeric(strsplit(main_log[j], split = " ")[[1]][length(strsplit(main_log[j], split = " ")[[1]])])
    }
    if(grepl("08 / 13", main_log[j])) {
      data$out_8[i] = as.numeric(strsplit(main_log[j], split = " ")[[1]][length(strsplit(main_log[j], split = " ")[[1]])])
    }
    if(grepl("09 / 13", main_log[j])) {
      data$out_9[i] = as.numeric(strsplit(main_log[j], split = " ")[[1]][length(strsplit(main_log[j], split = " ")[[1]])])
    }
    if(grepl("11 / 13", main_log[j])) {
      data$out_11[i] = as.numeric(strsplit(main_log[j], split = " ")[[1]][length(strsplit(main_log[j], split = " ")[[1]])])
    }
  }

  extra_log = readLines(file.path("log_files", log_fles[i]))
  for(j in 1 : length(extra_log)) {
    if(grepl("8 class start:", extra_log[j])) {
      data$f8_t0[i] = as.numeric(strsplit(extra_log[j], split = " ")[[1]][length(strsplit(extra_log[j], split = " ")[[1]])])
    }
    if(grepl("8 shift start:", extra_log[j])) {
      data$f8_t1[i] = as.numeric(strsplit(extra_log[j], split = " ")[[1]][length(strsplit(extra_log[j], split = " ")[[1]])])
    }
    if(grepl("8 finito:", extra_log[j])) {
      data$f8_tk[i] = as.numeric(strsplit(extra_log[j], split = " ")[[1]][length(strsplit(extra_log[j], split = " ")[[1]])])
    }
  }
  
}

write.csv(x = data, "assembly_times.csv")

setwd(dir = "./06_csvs/")
csvs_06 = list.files(full.names = TRUE)
setwd("../")

df_06 = data.frame(arr_size = vector(mode = "numeric"),
                     time1 = vector(mode = "numeric"),
                     time2 = vector(mode = "numeric"),
                     time3 = vector(mode = "numeric"),
                     time4 = vector(mode = "numeric"),
                     time5 = vector(mode = "numeric"),
                     time6 = vector(mode = "numeric"),
                     time7 = vector(mode = "numeric"),
                     time8 = vector(mode = "numeric"),
                     time9 = vector(mode = "numeric"),
                     time10 = vector(mode = "numeric"))
for(i in 1 : length(csvs_06)) {
  print(i)
  df_06 = rbind(df_06, t(read.csv(file = file.path("06_csvs", csvs_06[i]), header = T)))
}

write.csv(x = df_06, "csv06_times.csv")



setwd(dir = "./09_csvs/")
csvs_09 = list.files(full.names = TRUE)
setwd("../")

df_09 = data.frame(rep_no = vector(mode = "numeric"),
                   time1 = vector(mode = "numeric"),
                   time2 = vector(mode = "numeric"),
                   time3 = vector(mode = "numeric"),
                   time4 = vector(mode = "numeric"),
                   time5 = vector(mode = "numeric"),
                   time6 = vector(mode = "numeric"),
                   assembly = vector(mode = "character"),
                   i = vector(mode = "numeric"))
for(i in 1 : length(csvs_09)) {
  print(i)
  df_09[i, ] = c(t(read.csv(file = file.path("09_csvs", csvs_09[i]), header = T))[1,], "", 0)
  df_09$assembly[i] <- strsplit(csvs_09[i], split = "_")[[1]][length(strsplit(csvs_09[i], split = "_")[[1]]) - 3]
  df_09$i[i] <- strsplit(csvs_09[i], split = "_")[[1]][length(strsplit(csvs_09[i], split = "_")[[1]]) - 2]
}

write.csv(x = df_09, "csv09_times.csv")





