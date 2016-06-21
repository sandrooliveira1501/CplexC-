var generate = require('./generate_permutation')

for(var i = 4; i <= 8; i++){
  for (var j = 0; j < 5; j++){
    if(!generate.generatePermutation(i)){
      j--;
    }
  }
}
/*
var ord = generate.ordPermutation([0, 1, 2, 3, 4, 5, 9, 10, 11, 6, 7, 8, 15, 16, 17,12, 13, 14, 18, 19, 20, 21, 22]);

var ordOutput = "";

for(var i = 0; i < ord.length; i++){
  ordOutput += ord[i] + " ";
}

console.log(ordOutput);*/
