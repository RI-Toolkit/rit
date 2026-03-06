library(quantmod)

getSymbols("^AXJO", src = "yahoo", from = "1992-11-23",
           to = "2026-03-04")
asx200 <- Ad(AXJO)

df <- data.frame(date = index(asx200),
                 close = coredata(asx200))

df <- na.omit(df)

write.csv(df, "asx200_price.csv", row.names = FALSE)