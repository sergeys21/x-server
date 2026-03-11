// X0p_form fields enabling & disabling + fields validation
// Version-2.0 2013/05/19  fixed setOnloads as it did not work properly in IE after hitting "back" button

function LargeFloat(webpage, windowname, pagewidth, pageheight) {
  launch = window.open (webpage, windowname, 'width='+pagewidth+',height='+pageheight
                       +',toolbar=0,location=0,directories=0,status=0,menubar=0,scrollbars=1,resizable=1');
}

function setOnloads() {
    switchXrayFields();
}

function getXway() {
  for (var i=0; i < document.x0hform.xway.length; i++) {
    if (document.x0hform.xway[i].checked) {
      return document.x0hform.xway[i].value;
    }
  }
  return -1;
}

function switchXrayFields() {
  var xway = getXway();
  if (xway != 3) {
    document.x0hform.wave.disabled=false;
    document.x0hform.line.disabled=true;
  } else {
    document.x0hform.wave.disabled=true;
    document.x0hform.line.disabled=false;
  }
}

function x0h_validate() {
  var xway = getXway();
  if (xway == 1) {
    if (document.x0hform.wave.value == '') {
      alert('X-ray wavelength in Angstrem is not specified');
      return false;
    }
    var wave = parseFloat(document.x0hform.wave.value);
    if (wave <= 0.) {
      alert('X-ray wavelength='+wave+' Angstrom must be a positive value');
      return false;
    }
  }
  else if (xway == 2) {
    if (document.x0hform.wave.value == '') {
      alert('X-ray energy in KeV is not specified');
      return false;
    }
    var energy = parseFloat(document.x0hform.wave.value);
    if (energy <= 0.) {
      alert('X-ray energy='+energy+' KeV must be a positive value');
      return false;
    }
  }
  else if (xway == 3) {
    var selectedLine = document.x0hform.line.selectedIndex;
    if (document.x0hform.line.options[selectedLine].text == '') {
      alert('X-ray line selected, but not chosen from the list');
      return false;
    }
  }
  var selectedCode = document.x0hform.code.selectedIndex;
  if (document.x0hform.code.options[selectedCode].text == '') {
    alert('Crystal code is not selected');
    return false;
  }
// alert('Debug: OK'); return false;
  return true;
}
