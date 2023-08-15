TOPLOT <- plot1a

x <- c(2, 3, 4, 5, 10,25,50,100)
#x <- as.integer(rownames(plot2))
plot(x = x, y = TOPLOT$sire1, type = "b", ylim = c(0,1), xlim = c(2,100), log = "x",col = "red", xaxt='n',
     xlab = "Sample Size (log scale)", ylab = "Proportion containing MP", main = "Proportion of samples containing MP\nEqual sire proportions")
points(x = x, y = TOPLOT$sires2, type = "b", col = "blue")
points(x = x, y = TOPLOT$sires3, type = "b", col = "green")
points(x = x, y = TOPLOT$sires4, type = "b", col = "magenta")
legend(40, .5, legend=c("2 Sires", "3 Sires", "4 Sires"),
       col=c("blue","green","magenta"), lty = 1, cex=0.8)
axis(1, at=x)


TOPLOT <- plot1b

x <- c(2, 3, 4, 5, 10,25,50,100)
#x <- as.integer(rownames(plot2))
plot(x = x, y = TOPLOT$sire1, type = "b", ylim = c(0,1), xlim = c(2,100), log = "x",col = "red", xaxt='n',
     xlab = "Sample Size (log scale)", ylab = "Proportion containing MP", main = "Proportion of samples containing MP\nUnequal sire proportions")
points(x = x, y = TOPLOT$sires2, type = "b", col = "blue")
points(x = x, y = TOPLOT$sires3, type = "b", col = "green")
points(x = x, y = TOPLOT$sires4, type = "b", col = "magenta")
legend(40, .5, legend=c("2 Sires", "3 Sires", "4 Sires"),
       col=c("blue","green","magenta"), lty = 1, cex=0.8)
axis(1, at=x)





TOPLOT <- rbind(plot2a_ls,plot2a)

x <- c(2, 3, 4, 5, 10,25,50,100)
plot(x = x, y = TOPLOT$sires2, type = "b", ylim = c(0,1), col = "blue", log = "x", xaxt='n',
     xlab = "Sample Size (log scale)", ylab = "Proportion of iterations called positive for MP", 
     main = "Effect of Sample Size and Sire Number on Proportion of MP Calls\nEqual Sire, Equal Allele Proportions")
points(x = x, y = TOPLOT$sire1, type = "b", col = "red")
points(x = x, y = TOPLOT$sires3, type = "b", col = "green")
points(x = x, y = TOPLOT$sires4, type = "b", col = "magenta")
abline(h = .8, col = "gray")
legend(2, .95, legend=c("1 Sire","2 Sires", "3 Sires", "4 Sires"),
       col=c("red","blue","green","magenta"), lty = 1, cex=0.8)

TOPLOT <- rbind(plot2a_ls,plot2a)

x <- c(2, 3, 4, 5, 10,25,50,100)
plot(x = x, y = TOPLOT$sires2, type = "b", ylim = c(0,1), col = "blue", log = "x", xaxt='n',
     xlab = "Sample Size (log scale)", ylab = "Proportion of iterations called positive for MP", 
     main = "Effect of Sample Size and Sire Number on Proportion of MP Calls\nEqual Sire, Equal Allele Proportions")
points(x = x, y = TOPLOT$sire1, type = "b", col = "red")
points(x = x, y = TOPLOT$sires3, type = "b", col = "green")
points(x = x, y = TOPLOT$sires4, type = "b", col = "magenta")
abline(h = .8, col = "gray")
legend(2, .95, legend=c("1 Sire","2 Sires", "3 Sires", "4 Sires"),
       col=c("red","blue","green","magenta"), lty = 1, cex=0.8)
axis(1, at=x)

TOPLOT <- rbind(plot3b_smalls,plot3b)

x <- c(2, 3, 4, 5, 10,25,50,100)
plot(x = x, y = TOPLOT$sires2, type = "b", ylim = c(0,1), col = "blue", log = "x", xaxt='n',
     xlab = "Sample Size (log scale)", ylab = "Proportion of iterations called positive for MP", 
     main = "Effect of Sample Size and Sire Number on Proportion of MP Calls\nUnequal Sire, Equal Allele Proportions")
points(x = x, y = TOPLOT$sire1, type = "b", col = "red")
points(x = x, y = TOPLOT$sires3, type = "b", col = "green")
points(x = x, y = TOPLOT$sires4, type = "b", col = "magenta")
abline(h = .8, col = "gray")
legend(2, .95, legend=c("1 Sire","2 Sires", "3 Sires", "4 Sires"),
       col=c("red","blue","green","magenta"), lty = 1, cex=0.8)
axis(1, at=x)


TOPLOT <- rbind(plot3c_smalls,plot3c_bigs)

x <- c(2, 3, 4, 5, 10,25,50,100)
plot(x = x, y = TOPLOT$sires2, type = "b", ylim = c(0,1), col = "blue", log = "x", xaxt='n',
     xlab = "Sample Size (log scale)", ylab = "Proportion of iterations called positive for MP", 
     main = "Effect of Sample Size and Sire Number on Proportion of MP Calls\nEqual Sire, Unequal Allele Proportions")
points(x = x, y = TOPLOT$sire1, type = "b", col = "red")
points(x = x, y = TOPLOT$sires3, type = "b", col = "green")
points(x = x, y = TOPLOT$sires4, type = "b", col = "magenta")
abline(h = .8, col = "gray")
legend(2, .95, legend=c("1 Sire","2 Sires", "3 Sires", "4 Sires"),
       col=c("red","blue","green","magenta"), lty = 1, cex=0.8)
axis(1, at=x)


TOPLOT <- plot3d

x <- c(2, 3, 4, 5, 10,25,50,100)
plot(x = x, y = TOPLOT$sires2, type = "b", ylim = c(0,1), col = "blue", log = "x", xaxt='n',
     xlab = "Sample Size (log scale)", ylab = "Proportion of iterations called positive for MP", 
     main = "Effect of Sample Size and Sire Number on Proportion of MP Calls\nUnequal Sire, Unequal Allele Proportions")
points(x = x, y = TOPLOT$sire1, type = "b", col = "red")
points(x = x, y = TOPLOT$sires3, type = "b", col = "green")
points(x = x, y = TOPLOT$sires4, type = "b", col = "magenta")
abline(h = .8, col = "gray")
legend(2, .95, legend=c("1 Sire","2 Sires", "3 Sires", "4 Sires"),
       col=c("red","blue","green","magenta"), lty = 1, cex=0.8)
axis(1, at=x)




TOPLOT <- rbind(plot4_ls,plot4)

x <- c(1-.55, 1-.6, 1-.65, 1-.7, 1-.75, 1-.8, 1-.85, 1-.9, 1-.95, 1-.99)
plot(x = x, y = TOPLOT$sires2, type = "b", ylim = c(0,1), col = "blue", xaxt='n',
     xlab = "Target Type I error in full sib calls", ylab = "Proportion of iterations called positive for MP", 
     main = "Effect of Relationship Coefficient Cutoff on Proportion of MP Calls")
points(x = x, y = TOPLOT$sire1, type = "b", col = "red")
points(x = x, y = TOPLOT$sires3, type = "b", col = "green")
points(x = x, y = TOPLOT$sires4, type = "b", col = "magenta")
abline(h = .8, col = "gray")
legend(.35, .4, legend=c("1 Sire","2 Sires", "3 Sires", "4 Sires"),
       col=c("red","blue","green","magenta"), lty = 1, cex=0.8)
axis(1, at=x)



TOPLOT <- plot5

x <- c(2,3,4,5,6)
plot(x = x, y = TOPLOT$sires2, type = "b", ylim = c(0,1), col = "blue", xaxt='n',
     xlab = "Number of Alleles per Locus", ylab = "Proportion of iterations called positive for MP", 
     main = "Effect of Allelic Richness on Proportion of MP Calls")
points(x = x, y = TOPLOT$sire1, type = "b", col = "red")
points(x = x, y = TOPLOT$sires3, type = "b", col = "green")
points(x = x, y = TOPLOT$sires4, type = "b", col = "magenta")
abline(h = .8, col = "gray")
legend(5.5, .5, legend=c("1 Sire","2 Sires", "3 Sires", "4 Sires"),
       col=c("red","blue","green","magenta"), lty = 1, cex=0.8)
axis(1, at=x)


TOPLOT <- plot6

x <- c(5,10,15,20,25)
plot(x = x, y = TOPLOT$sires2, type = "b", ylim = c(0,1), col = "blue", xaxt='n',
     xlab = "Number of Loci", ylab = "Proportion of iterations called positive for MP", 
     main = "Effect of Locus Number on Proportion of MP Calls")
points(x = x, y = TOPLOT$sire1, type = "b", col = "red")
points(x = x, y = TOPLOT$sires3, type = "b", col = "green")
points(x = x, y = TOPLOT$sires4, type = "b", col = "magenta")
abline(h = .8, col = "gray")
legend(22, .6, legend=c("1 Sire","2 Sires", "3 Sires", "4 Sires"),
       col=c("red","blue","green","magenta"), lty = 1, cex=0.8)
axis(1, at=x)






