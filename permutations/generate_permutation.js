var _ = require('lodash');
var sh = require('shelljs');

var obj = {};
obj.generatePermutation = generatePermutation;

function generatePermutation(n){
  var permutation = [];
  var ord = [];

  while(permutation.length < n){

    var element =  Math.floor((Math.random() * n));

    if(!_.includes(permutation, element)){
      permutation.push(element);
    }

  }

  for(var i = 0; i < n; i++){

    ord[permutation[i]] = i;

  }

  //console.log(permutation);
  //console.log(ord);

  var entrada = "" + (permutation[0] + 1);
  for(var i = 1; i < permutation.length; i++){
    entrada += "," + (permutation[i] + 1);
  }
  //console.log(entrada);
  var command = "python2.6 ../code_dias2012/elias2005.py " + entrada;
  //console.log(command);
  var process = sh.exec(command, {silent:true});
  var output = _.split(process.stdout, " - ");
  var l = output[2];
  //console.log(l);

  var permutationOutput = "";
  var ordOutput = "";

  for(var i = 0; i < n; i++){
    permutationOutput += permutation[i] + " ";
    ordOutput += ord[i] + " ";
  }

  permutationOutput = _.trimEnd(permutationOutput);
  ordOutput = _.trimEnd(ordOutput);

  console.log(permutationOutput)
  console.log(ordOutput)
  console.log(l)
}

module.exports = obj;
