
# Get function help. In R Studio you can also press F1 button when cursor is in the function.
?data()

data()


iris
?iris



# Inspecting data ---------------------------------------------------------

# Show first six rows
head(iris)
head(iris, 3)
# Show last six rows
tail(iris)
# Show object's dimension
dim(iris)
# Show number of rows
nrow(iris)
# Object class
class(iris)
# Object structure
str(iris)

View(iris)
# Subsetting 
iris[1, 1]
iris[1, 1] <- 2

# Data types --------------------------------------------------------------

## Numeric
# Printing one number or sequence of numbers
1
1:10
# Combining values
c(1, 5)
# Create number sequence and store into variable
numbers <- seq(0, 10, 2.5)
numbers
# Mathematical operations
numbers + 1
numbers * 1
# Statistical operations
mean(numbers)
sum(numbers)
sd(numbers)
# Check object class
class(numbers)
# In R float and integer are both numeric
class(1)
class(1.0)

## String
strings <- "numbers"

class(strings)

library(dplyr)

## Boolean
# Less than
3 < 5
# Equals
3 == 3
# Not equals
3 != 3

## Logical with boolean
3 & 5 < 8
3 & 10 < 8

"a" < "b"
"c" < "b"

# Special
# Doesn't exist
?NULL
# Missing value
?NA
# Not a number
?NaN

# Inspecting data ---------------------------------------------------------
## This part might be merged with the above
# Clear environment 
rm(list = ls())

install.packages("dplyr")
install.packages("magrittr")

library(dplyr)
library(magrittr)

iris
str(iris)
summary(iris)

table(iris$Species)
unique(iris$Species)
length(iris$Species)

iris$Species

iris$Species == "setosa"
iris[iris$Species == "setosa", ]

subset(iris, Species == "setosa")
iris_setosa <- subset(iris, Species == "setosa")
iris_setosa_large_sepal <- iris_setosa[iris_setosa$Sepal.Length > 5, ]
nrow(iris_setosa_large_sepal)

# With dplyr

# %>% 
iris %>% 