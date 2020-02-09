# https://www.nltk.org/book/ch01.html
import string
import nltk
from nltk.tokenize import word_tokenize

nltk.download('book')
nltk.download('punkt')

from nltk.book import text1, text2, text3, text4, text5, text6, text7, text8, text9

text1
text1.concordance('monstrous')
len(text1)

comment = 'With the Senate trial of Mr. Trump now underway, we deployed a team of journalists to find out. We contacted hundreds of voters who had responded to an online survey saying they would be willing to be interviewed. We reached 81 people, from nearly 30 states. They were Democrats, Republicans and independents. They were retirees and real estate agents, teachers and stay-at-home parents. The youngest was 21; the oldest was 82. Even before the opening statements at the trial had begun, most had already made up their minds on their preferred verdict. As one independent voter from Ohio put it, “maybe they should ask the people what they should do. It should be our vote.”'

def count_words(s):
    s = str(s).lower()
    s = s.translate(str.maketrans('', '', string.punctuation))
    return len(word_tokenize(s))

count_words(comment)

# Average number of words
