# Compute the independence model for the pewee data.
Indep <- march.indep.construct(pewee)
# Display the model
print(Indep)
# Compute the half-length 95% confidence interval for each element of the distribution.
march.indep.bailey(Indep,alpha=0.05)

# Compute a second-order MTDg model for the pewee data.
MTD2g <- march.mtd.construct(pewee,2,mtdg=TRUE)
# Display the model
print(MTD2g)
# Compute the half-length 95% confidence interval for all parameters
# of the MTD2g model.
march.mtd.bailey(MTD2g,alpha=0.05)
