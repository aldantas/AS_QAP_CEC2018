<table>
  <tbody>
    <tr>
      <th>File</th>
      <th align="center">Description</th>
    </tr>
    <tr>
      <td><code>performance.log</code></td>
      <td align="left">Textual output reporting the overall Accuracy and the
      individual classes performances (Precision, Recall and F-Score)</td>
    </tr>
    <tr>
      <td><code>predicts.dat</code></td>
      <td align="left">A <a
      href="https://docs.python.org/3/library/pickle.html">pickle</a>
      file storing a tuple with five uni-demensional lists, each containing
      information of all 10 folds, which are:
	 <ul>
	  <li>True labels</li>
	  <li>Predicted labels</li>
	  <li>Probability estimates</li>
	  <li>Scores by fold</li>
	  <li>Features importance by fold</li>
	</ul>
      </td>
    </tr>
  </tbody>
</table>
