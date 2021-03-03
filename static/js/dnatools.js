/**
 * Fetches some data from some URL given some input
 * @param i a number
 * @param selected another number
 */
let get_new_plot = function(i, selected) {


    let url = new URL("/get_details", location.href);
    let params = {
            source: i,
            index: selected
        };
    // Combine URL and parameters to format a GET-request URL (i.e.: http://localhost:8080/get_details?source=2&index=4)
    Object.keys(params).forEach(key => url.searchParams.append(key, params[key]));

    fetch(url)
        .then((resp) => resp.json())
        .then(function(data) {
            // Print what we've received. Update this to receive a new plot object and use 'Bokeh.embed' to include it.
            console.log("Data received:");
            console.log(data);
        })
};

$(function () {

    /**
     * Javascript functions and events
     */

    $("#gc-submit").click(function(event) {
        // Make sure the form doesn't actually submit
        event.preventDefault();

        // Get the selected file
        const formData = new FormData();
        const fileField = document.querySelector('input[type="file"]');
        formData.append('file', fileField.files[0]);

        // TODO: check if a file is selected, ignore if not
        
        // Get the plot data given the uploaded file
        fetch('/gc-percentage', {
          method: 'POST',
          body: formData,
        })
          .then(function (response) {
            return response.json();
          })
          .then(function (plot) {
            // Empty the div (otherwise plots are appended)
            $("#GCplot").html("");
            // Embed the plot in the given div
            Bokeh.embed.embed_item(plot, "GCplot");
          });
      });
});