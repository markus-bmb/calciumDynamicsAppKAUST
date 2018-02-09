apTrain = ActionPotentialTrain(0.0, 1.0, 50, -65.0)

t = 0
while t < 0.2 do
	print(t .. "  " .. apTrain:membrane_potential(t))
	t = t + 1e-5
end

