## script does not exist yet

## fragments here were just pulled out of the main DEA script

data.wide <- cbind(expDes, t(y$counts[rownames(tt)[1:4], ]))
data.tall <- melt(data.wide,
                  id.vars = names(expDes),
                  variable.name = 'contig', value.name = 'Expression')
str(data.tall)

#x <- data.frame(expDes, gExp = y$counts[rownames(tt)[1], ])
p <- ggplot(data.tall, aes(x = gType, y = Expression))
p + geom_point() + facet_grid(tx ~ contig)
