import numpy as np

def bootstrap (func, data, sign, weight, maxNumBootstrapIteration):
  ####
  #### convert to numpy structures
  ####
  data   = np.array(data)
  sign   = np.array(sign)
  weight = np.array(weight)

  # combine weights and sign
  signAndWeight = sign * weight

  # compute true mean
  trueDataMean = np.tensordot(data, signAndWeight, axes=([1, 0])) / np.sum(signAndWeight)
  trueMean     = func(*trueDataMean)

  # bootstrap error
  means = []
  for j in range(0, min(weight.size - 1, maxNumBootstrapIteration)):
    # pick data points on random
    I = np.random.randint(0, weight.size - 1, weight.size)
    _data          = data[:, I]
    _signAndWeight = signAndWeight[I]

    # compute data mean value
    dataMean = np.tensordot(_data, _signAndWeight, axes=([1, 0])) / np.sum(_signAndWeight)

    # function of data mean value
    means += [func(*dataMean)]

  # estimate error
  if len(means): err = np.mean((means - trueMean)**2)**0.5
  else:          err = 0

  return [trueMean, err]