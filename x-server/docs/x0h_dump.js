
function LargeFloat(webpage, windowname, pagewidth, pageheight) {
  launch = window.open (webpage, windowname, 'width='+pagewidth+',height='+pageheight
                       +',toolbar=0,location=0,directories=0,status=0,menubar=0,scrollbars=1,resizable=1');
}

function x0h_validate() {
    var selectedCode = document.x0hdump.code.selectedIndex;
    if (document.x0hdump.code.options[selectedCode].text == '') {
      alert('Material code is not selected');
      return false;
    }
// alert('Debug: OK'); return false;
   return true;
}
