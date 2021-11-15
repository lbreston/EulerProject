# Signal processing library 
using DSP

# Generating function for distribution of sums:  f(s) = ΣG(x,y) = Σg_n,k*x^n*y^k = ∏_n(1+x*y^k) 
# where n is the subset size and k is value of the sum 

# creates vector of factors of generator function P[k] = (1+x*y^k) for sequence of squares from 1:k_max²
function generating_factors(k_max)
    #* n+1 x k+1 matrix encoding: row = size of subset, column = degree of value of sum
    #* first row and column represent n,k=0 
    function f̨ₖ(k)
        f=zeros(Int8,2,k^2+1)
        f[1,1]=1
        f[2,k^2+1]=1
        return f
    end
    return [f̨ₖ(k) for k in 1:k_max]
end

# computes coefficients g_n_k for subsets of size n
function compute_generating_coefs(G, n)
    function clipvals(x)
        x=abs(round(x))
        if x==0
            return 0
        elseif x==1
            return 1
        else
            return 2
        end
    end
   # computes product using convolutions. Convolutions use FFTs to make it go way faster!
    coefs = reduce((x,y)->clipvals.(conv(x,y)), G)
    vals =  collect(0:size(coefs,2)-1)
    return coefs[n+1, :], vals 
end

# sum of the unique subset sums 
function sum_unique_subset_sums(k_max, n)
    G = generating_factors(k_max)
    coefs, vals = compute_generating_coefs(G, n)
    return sum(vals[coefs.==1])
end


sum_unique_subset_sums(100,50)



