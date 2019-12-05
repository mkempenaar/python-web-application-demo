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
        formData.append('name', 'Henk');
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