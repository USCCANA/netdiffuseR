# Set dimensions
n <- 4  # Number of nodes
t <- 3  # Number of time steps
q <- 2  # Number of contagions

# Define graph (n x n matrix)
graph <- matrix(c(0, 1, 0, 0,
                  1, 0, 1, 0,
                  0, 1, 0, 1,
                  0, 0, 1, 0), nrow = n, byrow = TRUE)

# Define attrs (n x t matrix)
attrs <- matrix(c(1, 2, 3,
                  4, 5, 6,
                  7, 8, 9,
                  10,11,12), nrow = n)
#attrs <- array(c(1,2,3,
#                 4,5,6,
#                 7,8,9,
#                 10,11,12,
# Second contagion
#                 2,4,6,
#                 8,10,12,
#                 14,16,18,
#                 20,22,24), dim = c(n,t,q))

# Define cumadopt (n x t x q array)
cumadopt <- array(c(1,2,3,
                    4,5,6,
                    7,8,9,
                    10,11,12), dim = c(n,t))

cumadopt <- array(c(1,2,3,
                    4,5,6,
                    7,8,9,
                    10,11,12,
                    # Second contagion
                    2,4,6,
                    8,10,12,
                    14,16,18,
                    20,22,24), dim = c(n,t,q))

# ORIGINAL: Graph -> n x n
#           attrs -> n x T
#           cumadopt-> n x T
#           ans   -> 1 x n
#
# so (attrs * cumadopt) -> n x T
#    (graph %*% (attrs * cumadopt)) -> n x T
#
# normalization: as.vector(ans/( graph %*% attrs + 1e-20 ))
# so (graph %*% attrs) -> n x T

# NEW:      Graph -> n x n
#           attrs -> n x T
#           cumadopt-> n x T x q
#           ans   -> q x n
#
# so (graph %*% (attrs * cumadopt)) -> n x T x q
#
# normalization
# so (graph %*% attrs) -> n x T

# ANS

ans <- apply(cumadopt, MARGIN=3, function(ca) graph %*% (attrs * ca))
ans <- array(ans, dim = c(n,t,q))
dim(ans)

# NORMALIZATION

den <- graph %*% attrs
dim(den)

ans/( graph %*% attrs + 1e-20 ) # Error

# TRY THIS

ans_norm <- apply(cumadopt, MARGIN=3, function(ca) graph %*% (attrs * ca) / ( graph %*% attrs + 1e-20 ))
ans <- array(ans, dim = c(n,t,q))
dim(ans)

# OR THIS

#initializing array
ans <- array(0, dim = c(n,t,q))
norm <- graph %*% attrs + 1e-20

normalized <- TRUE

#loop for q contagions
for (k in seq_len(q)) {
  if (normalized) ans[,,k] <- graph %*% (attrs * cumadopt[,,k]) / norm
  else ans[,,k] <- graph %*% (attrs * cumadopt[,,k])
}

as.vector(ans, dim = c(n,t,q))

dim(ans)


##

ans <- ( graph %*% (attrs * cumadopt) )
dim(ans) # n x t

ans_norm <- ans/( graph %*% attrs + 1e-20 )
class(ans_norm) # "matrix" "array"
dim(ans_norm) # n x t

ans_norm_vec <- as.vector(ans/( graph %*% attrs + 1e-20 ))
class(ans_norm_vec) # "numeric"
dim(ans_norm_vec) # NULL, only a vector o length 4x3=12
