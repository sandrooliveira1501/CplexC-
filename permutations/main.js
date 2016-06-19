var generate = require('./generate_permutation')

for(var i = 8; i <= 12; i++){
  for (var j = 0; j < 5; j++){
    if(!generate.generatePermutation(i)){
      j--;
    }
  }
}
