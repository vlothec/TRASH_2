# https://stackoverflow.com/questions/27312292/convert-seconds-to-days-hoursminutesseconds @the_10th_letter_x2
dhms <- function(t) {
  paste(t %/% (60 * 60 * 24),
        paste(formatC(t %/% (60 * 60) %% 24, width = 2, format = "d", flag = "0"),
              formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0"),
              formatC(t %% 60, width = 2, format = "d", flag = "0"),
              sep = ":"))
}