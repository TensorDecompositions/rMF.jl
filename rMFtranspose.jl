"Transpose non-numeric vector"
function transposevector(a)
	reshape(a, 1, length(a))
end

"Transpose non-numeric matrix"
function transposematrix(a)
	permutedims(a, (2, 1))
end