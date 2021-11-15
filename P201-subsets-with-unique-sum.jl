using DSP

#* Generating function for distribution of sums:  f(s) = ΣG(x,y) = Σg_n,k*x^n*y^k = ∏_n(1+x*y^k) 
#* where n is the subset size and k is value of the sum 

#* n+1 x k+1 matrix encoding: row = size of subset, column = degree of value of sum
#* first row and column represent n,k=0 
function pfun(d,n)
    z=zeros(Float64,2,d^2+1)
    z[1,1]=1
    z[2,n^2+1]=1
    return z
end

#* vector of factors of generator function P[k] = (1+x*y^k) 
degree=100
P=[pfun(degree,k) for k in 1:degree]

#* compute product using convolutions. Convolutions use FFTs to make it go way faster!
G=reduce((x,y)->abs.(round.(conv(x,y))), P)

#* Find sums of 50 element subsets that appear once i.e. g_51,k = 1 
idx=findall(G[51,:].==1)

#* Sum unique sums 
sumval=0:1000001
sum(sumval[idx])
