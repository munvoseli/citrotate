
--[[
trying this:
https://en.wikipedia.org/wiki/Cross_product#Coordinate_notation
but, say, cross(ijk)=l instead of i cross j = k
and also applying anticommutativity with swaps

oops, no longer trying that, thanks quintopia
--]]


function dot(a, b)
	local sum = 0
	for i=1,#a do
		sum = sum + a[i] * b[i]
	end
	return sum
end
function scmul(x, v)
	local v2 = {}
	for i=1,#v do
		table.insert(v2, x * v[i])
	end
	return v2
end
function proj(a, u)
	return scmul(dot(a, u) / dot(u, u), u)
end
function vsub(a, b)
	local v2 = {}
	for i=1,#a do
		table.insert(v2, a[i] - b[i])
	end
	return v2
end
function normeq(v) -- scales magnitude to 1
	local mag = 0
	for i=1,#v do
		mag = mag + v[i] * v[i]
	end
	mag = math.sqrt(mag)
	for i=1,#v do
		v[i] = v[i] / mag
	end
end

function gramschmidt(vecs, n, v) -- 1 <= n < #vecs
	local v = scmul(1, v) -- just for copy
	normeq(v)
	for i=1,n do
		v = vsub(v, proj(v, vecs[i]))
		normeq(v)
	end
	return v
end

-- for determining which axis unit vector
-- is most similar to the given vector
function maxcomponent(v)
	local x = 0
	local k = 0
	for i=1,#v do
		if math.abs(v[i]) > x then
			x = math.abs(v[i])
			k = i
		end
	end
	return k
end
function getcompvec(c, n)
	local v = {}
	for i=1,n do
		table.insert(v, 0)
	end
	v[c] = 1
	return v
end

function orthonormal(vec)
	local vec = scmul(1, vec)
	normeq(vec)
	local vecs = {vec} -- column vectors
	local c = maxcomponent(vec)
	--[[
	idea is that, it should, before gram schmidt,
	be roughly equivalent to identity matrix.
	if the provided vector is kinda negative unit, that flips.
	also, flips based on position of vec (c):
		c == 1, vec in front: 0 swaps -> 0 flips
		c == 2: 1 swap -> flip
		c == 3: 2 swaps -> noflip
		c even -> flip
		c  odd -> noflip
	thus, guessing that flip == vec[c]<0 xor c%2==0
	--]]
	function xor(a, b)
		return not(a == b)
	end
	local shouldFlip = xor(vec[c] < 0, c % 2 == 0)
	for i=1,c-1 do
		table.insert(vecs, getcompvec(i, #vec))
	end
	for i=c+1,#vec do
		table.insert(vecs, getcompvec(i, #vec))
	end
	if shouldFlip then
		for i=1,#vec do
			vecs[2][i] = vecs[2][i] * -1
		end
	end
	--printSquareMatrix(vecs)
	for i=2,#vecs do
		vecs[i] = gramschmidt(vecs, i-1, vecs[i])
	end
	return vecs
end

function printSquareMatrix(colvecs)
	for row=1,#colvecs do
		str = ""
		for col=1,#colvecs do
			str = str .. colvecs[col][row] .. "\t"
		end
		str = str .. "\n"
		print(str)
	end
end


function testMatrix(colvecs)
	-- check that each column is orthogonal to every other
	print("Should all be nearly 0:")
	for i=1,#colvecs do
	for j=1,#colvecs do
	if j ~= i then
		print(dot(colvecs[i], colvecs[j]))
	end
	end
	end
	print("Should all be nearly 1:")
	for i=1,#colvecs do
		print(dot(colvecs[i], colvecs[i]))
	end
end

--local mat = orthonormal({1, 2, 0, 0, 0})
--printSquareMatrix(mat)
--testMatrix(mat)

local M = {}

--[[
accepts vector
(internally, will scale that to length 1 if it's not already 1).
returns array of column vectors
  which represent orthonormal (rotation) matrix.

or, hopefully it's a rotation matrix.
it might flip things sometimes.
orthonormal matrices can flip things.
if it does that, please contact me,
i have no idea how to fix it atm.
--]]
M.getOrthonormalMatrix = orthonormal

return M
