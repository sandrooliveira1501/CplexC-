var _ = require('lodash');
var sh = require('shelljs');

var obj = {};
obj.generatePermutation = generatePermutation;
obj.ordPermutation = ordPermutation;

function ordPermutation(permutation){
  var ord = [];
  for(var i = 0; i < permutation.length; i++){

    ord[permutation[i]] = i;

  }

  return ord;
}

function generatePermutation(n, rep){
  var permutation = [];
  var ord = [];

  while(permutation.length < n){

    var element =  Math.floor((Math.random() * n));

    if(!_.includes(permutation, element)){
      permutation.push(element);
    }

  }

  var newN = n * rep;
  var newPermutation = [];
  for(var i = 0; i < n; i++){
    var elem = permutation[i] * rep;
    for(var e = 0; e < rep; e++){
      newPermutation.push(elem + e);
    }
  }

  permutation = newPermutation;
  n = newN;

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
  var command = "python2.6 ../code_dias2012/dias2010_1_375.py " + entrada;
  var command2 = "python2.6 ../code_dias2012/elias2005.py " + entrada;
  var command3 = "python2.6 ../code_dias2012/dias2010_1_5.py " + entrada;

  var process = sh.exec(command, {silent:true});
  var output = _.split(process.stdout, " - ");
  var l = output[2];

  if(output[0].startsWith("Erro")){
    process = sh.exec(command2, {silent:true});
    output = _.split(process.stdout, " - ");
    l = output[2];

    if(output[0].startsWith("Erro")){
      process = sh.exec(command3, {silent:true});
      output = _.split(process.stdout, " - ");
      l = output[2];

      if(output[0].startsWith("Erro")){
        return false;
      }
    }

  }

  var permutationOutput = "";
  var ordOutput = "";

  for(var i = 0; i < n; i++){
    permutationOutput += permutation[i] + " ";
    ordOutput += ord[i] + " ";
  }

  permutationOutput = _.trimEnd(permutationOutput);
  ordOutput = _.trimEnd(ordOutput);

  console.log(permutationOutput);
  console.log(ordOutput);
  if(l == 0){
    l = 1;
  }
  console.log(parseInt(l));

  return true;
}

module.exports = obj;
