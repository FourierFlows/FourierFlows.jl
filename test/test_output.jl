function testwithoutjld2()
  namewithjld2 = "blahblah.jld2"
  namewithoutjld2 = "blahblah"
  
  FourierFlows.withoutjld2(namewithjld2) == namewithoutjld2 && FourierFlows.withoutjld2(namewithoutjld2) == namewithoutjld2
end
