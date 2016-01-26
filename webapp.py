from flask import Flask, render_template, request
from logic import translate_dna

app = Flask(__name__)


@app.route('/')
def hello_world():
    """ Renders the DNA tools web form
    :return: the rendered view
    """
    return render_template('dnatools.html', title='Translate DNA Sequence')


@app.route('/dna-translate', methods=['POST'])
def process_webform():
    """ Receives a DNA sequence and selected frames and performs translation.
    The resulting dictionary with 'frame:translated sequence' is passed to the
    template and rendered as a table.
    :return: the rendered view containing the result table
    """
    # Retrieve the entered DNA-sequence (<textarea>)
    sequence = request.form.get('seq')

    # Retrieve the requested frames (<checkbox>; both forward and reverse)
    # and cast each list element to integer values
    fframes = [int(frame) for frame in request.form.getlist('fframe')]
    rframes = [int(frame) for frame in request.form.getlist('rframe')]

    # Translate the sequence with the requested frames
    translated_forward_seqs = translate_dna(sequence, fframes)
    translated_reverse_seqs = translate_dna(sequence, rframes, reverse=True)

    # Combine dictionaries
    translated_seqs = {**translated_forward_seqs, **translated_reverse_seqs}

    # Render result page
    return render_template('translated.html',
                           title="Translated DNA Sequence",
                           translated=translated_seqs)

if __name__ == '__main__':
    app.run()
